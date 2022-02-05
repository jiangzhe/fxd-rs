use crate::utils::get_units;
use std::hash::{Hash, Hasher};
use crate::error::{Result, Error};

pub const MAX_UNITS: usize = 9;
pub const MAX_FRAC: usize = 30;
pub const MAX_DIGITS: usize = 65;
pub const DIGITS_PER_UNIT: usize = 9;
pub const UNIT: i32 = 1_000_000_000;
pub const HALF_UNIT: i32 = UNIT / 2;
pub const MAX_FRAC_UNITS: usize = (MAX_FRAC + DIGITS_PER_UNIT - 1) / DIGITS_PER_UNIT;
pub const MAX_I64: FixedDecimal = FixedDecimal {
    intg: 19 as i8,
    frac: 0,
    lsu: [854775807, 223372036, 9, 0, 0, 0, 0, 0, 0],
};
pub const MIN_I64: FixedDecimal = FixedDecimal {
    intg: (19 | 0x80) as i8,
    frac: 0,
    lsu: [854775808, 223372036, 9, 0, 0, 0, 0, 0, 0],
};
pub const MAX_U64: FixedDecimal = FixedDecimal {
    intg: 20 as i8,
    frac: 0,
    lsu: [709551615, 446744073, 18, 0, 0, 0, 0, 0, 0],
};

/// FixedDecimal is an implementation of "exact number" type defined
/// SQL standard.
/// It directly maps to MySQL data type "DECIMAL".
/// it has fixed allocation: 9*4+2=38 bytes.
#[derive(Debug, Clone)]
pub struct FixedDecimal {
    // Integral digits.
    // Maximum is 65, if frac is 0.
    // Use most significant bit to indicates whether it is negative
    pub(crate) intg: i8,
    // Fractional digits.
    // Maximum is 30.
    pub(crate) frac: i8,
    pub(crate) lsu: [i32; MAX_UNITS],
}

impl FixedDecimal {
    #[inline]
    pub fn zero() -> FixedDecimal {
        FixedDecimal {
            intg: 1,
            frac: 0,
            lsu: [0; MAX_UNITS],
        }
    }

    #[inline]
    pub fn one() -> FixedDecimal {
        let mut fd = Self::zero();
        fd.lsu[0] = 1;
        fd
    }

    #[inline]
    pub fn is_neg(&self) -> bool {
        self.intg & (0x80u8 as i8) != 0
    }

    #[inline]
    pub fn set_neg(&mut self) {
        self.intg |= 0x80u8 as i8;
    }

    #[inline]
    pub(crate) fn set_neg_and_check_zero(&mut self) {
        self.set_neg();
        if self.is_zero() {
            self.set_pos();
        }
    }

    #[inline]
    pub fn set_pos(&mut self) {
        self.intg &= 0x7f;
    }

    #[inline]
    pub fn set_zero(&mut self) {
        self.intg = 1;
        self.frac = 0;
        self.reset_units();
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.lsu.iter().all(|n| *n == 0)
    }

    #[inline]
    pub fn set_one(&mut self) {
        self.intg = 1;
        self.frac = 0;
        self.reset_units();
        self.lsu[0] = 1;
    }

    #[inline]
    pub fn intg(&self) -> i8 {
        self.intg & 0x7f
    }

    #[inline]
    pub fn intg_units(&self) -> usize {
        get_units(self.intg())
    }

    #[inline]
    pub fn frac(&self) -> i8 {
        self.frac
    }

    #[inline]
    pub fn frac_units(&self) -> usize {
        get_units(self.frac())
    }

    #[inline]
    fn reset_units(&mut self) {
        self.lsu.iter_mut().for_each(|n| *n = 0);
    }

    #[inline]
    pub fn as_i64(&self) -> Result<i64> {
        if self < &MIN_I64 || self > &MAX_I64 {
            return Err(Error::ExceedsConversionTargetRange)
        }
        if self.frac_units() > 0 {
            let mut target = Self::zero();
            self.round_to(&mut target, 0);
            return Ok(target.lsu_to_i64())
        }
        Ok(self.lsu_to_i64())
    }

    fn lsu_to_i64(&self) -> i64 {
        if self.is_neg() {
            -self.lsu[0] as i64 - self.lsu[1] as i64 * UNIT as i64 - self.lsu[2] as i64 * UNIT as i64 * UNIT as i64
        } else {
            self.lsu[0] as i64 + self.lsu[1] as i64 * UNIT as i64 + self.lsu[2] as i64 * UNIT as i64 * UNIT as i64
        }
    }

    #[inline]
    pub fn as_u64(&self) -> Result<u64> {
        if self.is_neg() || self > &MAX_U64 {
            return Err(Error::ExceedsConversionTargetRange)
        }
        if self.frac_units() > 0 {
            let mut target = Self::zero();
            self.round_to(&mut target, 0);
            return Ok(target.lsu_to_u64())
        }
        Ok(self.lsu_to_u64())
    }

    fn lsu_to_u64(&self) -> u64 {
        self.lsu[0] as u64 + self.lsu[1] as u64 * UNIT as u64 + self.lsu[2] as u64 * UNIT as u64 * UNIT as u64
    }
}

/// Simple hash implementation on FixedDecimal.
///
/// NOTE: Same numeric values with different number of
/// internal units will have different hash code.
impl Hash for FixedDecimal {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // hash neg flag
        state.write_u8(if self.is_neg() { 1 } else { 0 });
        // hash nbr of intg units
        let intg_units = self.intg_units();
        state.write_u8(intg_units as u8);
        // hash nbr of frac units
        let frac_units = self.frac_units();
        state.write_u8(frac_units as u8);
        // hash all units
        for v in &self.lsu[..intg_units + frac_units] {
            state.write_u32(*v as u32)
        }
    }
}

impl From<i64> for FixedDecimal {
    fn from(src: i64) -> Self {
        if src == i64::MIN {
            return MIN_I64.clone();
        }
        if src < 0 {
            let mut fd = FixedDecimal::from(-src as u64);
            fd.set_neg();
            fd
        } else {
            FixedDecimal::from(src as u64)
        }
    }
}

impl From<u64> for FixedDecimal {
    fn from(mut src: u64) -> Self {
        let mut fd = FixedDecimal::zero();
        let mut i = 0;
        while src != 0 {
            let q = src / UNIT as u64;
            let r = src - q * UNIT as u64;
            fd.lsu[i] = r as i32;
            i += 1;
            src = q;
        }
        fd.intg = (i * DIGITS_PER_UNIT) as i8;
        fd
    }
}

#[cfg(test)]
mod tests {

    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_one() {
        let fd = FixedDecimal::one();
        assert!(!fd.is_neg());
        assert!(!fd.is_zero());
    }

    #[test]
    fn test_from_str_rand() {
        use rand::prelude::*;
        let mut rng = rand::thread_rng();
        let mut fd = FixedDecimal::zero();
        for _ in 0..128 {
            let intg: u32 = rng.gen_range(0..1 << 30);
            let frac: u32 = rng.gen_range(0..1 << 20);
            let s = if frac > 0 {
                format!("{}.{}", intg, frac)
            } else {
                format!("{}", intg)
            };
            assert!(fd.from_ascii_str(&s, true).is_ok());
            println!("s={}, fd={:?}", s, fd);
        }
    }

    #[test]
    fn test_hash() {
        for (s1, s2) in vec![("0.1", "0.10"), ("1.2", "1.2")] {
            let fd1: FixedDecimal = s1.parse().unwrap();
            let fd2: FixedDecimal = s2.parse().unwrap();
            assert_eq!(decimal_hash(&fd1), decimal_hash(&fd2));
        }
    }

    #[test]
    fn test_from_i64() {
        for (i, s) in vec![
            (0i64, "0"),
            (1, "1"),
            (-1, "-1"),
            (100, "100"),
            (-100, "-100"),
            (12345123456789, "12345123456789"),
            (-12345123456789, "-12345123456789"),
            (-9223372036854775808, "-9223372036854775808"),
            (9223372036854775807, "9223372036854775807"),
        ] {
            let fd1 = FixedDecimal::from(i);
            let fd2 = FixedDecimal::from_str(s).unwrap();
            assert_eq!(fd1, fd2);
            println!("{:?}", fd1);
        }
    }

    #[test]
    fn test_to_i64() {
        // success
        for (s, i) in vec![
            ("0", 0i64),
            ("1", 1),
            ("1.1", 1),
            ("1.5", 2),
            ("3.8", 4),
            ("-1.2", -1),
            ("-5.9", -6),
            ("100", 100),
            ("12345123456789", 12345123456789),
            ("-9223372036854775808", -9223372036854775808),
            ("9223372036854775807", 9223372036854775807),
        ] {
            let fd = FixedDecimal::from_str(s).unwrap();
            let res = fd.as_i64().unwrap();
            assert_eq!(res, i)
        }
        // fail
        for s in vec![
            "-9223372036854775809",
            "9223372036854775808",
            "10000000000000000000000",
        ] {
            let fd = FixedDecimal::from_str(s).unwrap();
            assert!(fd.as_i64().is_err());
        }
    }

    #[test]
    fn test_from_u64() {
        for (i, s) in vec![
            (0u64, "0"),
            (1, "1"),
            (100, "100"),
            (12345123456789, "12345123456789"),
            (18446744073709551615, "18446744073709551615"),
        ] {
            let fd1 = FixedDecimal::from(i);
            let fd2 = FixedDecimal::from_str(s).unwrap();
            assert_eq!(fd1, fd2);
            println!("{:?}", fd1);
        }
    }

    #[test]
    fn test_to_u64() {
        // success
        for (s, i) in vec![
            ("0", 0u64),
            ("1", 1),
            ("1.1", 1),
            ("1.5", 2),
            ("3.8", 4),
            ("100", 100),
            ("12345123456789", 12345123456789),
            ("9223372036854775807", 9223372036854775807),
            ("18446744073709551615", 18446744073709551615),
        ] {
            let fd = FixedDecimal::from_str(s).unwrap();
            let res = fd.as_u64().unwrap();
            assert_eq!(res, i)
        }
        // fail
        for s in vec![
            "-1",
            "18446744073709551616",
            "10000000000000000000000",
        ] {
            let fd = FixedDecimal::from_str(s).unwrap();
            assert!(fd.as_u64().is_err());
        }
    }

    #[inline]
    fn decimal_hash(fd: &FixedDecimal) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        let mut hasher = DefaultHasher::new();
        fd.hash(&mut hasher);
        hasher.finish()
    }
}
