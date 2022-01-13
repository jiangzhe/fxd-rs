use crate::utils::get_units;

pub const MAX_UNITS: usize = 9;
pub const MAX_FRAC: usize = 30;
pub const MAX_DIGITS: usize = 65;
pub const DIGITS_PER_UNIT: usize = 9;
pub const UNIT: i32 = 1_000_000_000;
pub const HALF_UNIT: i32 = UNIT / 2;
pub const MAX_FRAC_UNITS: usize = (MAX_FRAC + DIGITS_PER_UNIT - 1) / DIGITS_PER_UNIT;


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
        FixedDecimal{intg: 1, frac: 0, lsu: [0; MAX_UNITS]}
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
        self.frac= 0;
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
}


#[cfg(test)]
mod tests {

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
            let intg: u32 = rng.gen_range(0..1<<30);
            let frac: u32 = rng.gen_range(0..1<<20);
            let s = if frac > 0 {
                format!("{}.{}", intg, frac)
            } else {
                format!("{}", intg)
            };
            assert!(fd.from_ascii_str(&s, true).is_ok());
            println!("s={}, fd={:?}", s, fd);
        }
    }
}