use crate::arith::add_with_carry;
use crate::decimal::{FixedDecimal, DIGITS_PER_UNIT, HALF_UNIT, UNIT};
use crate::utils::{div9, get_units, mod9};

impl FixedDecimal {
    /// Round rounds this decimal with provided frac.
    /// frac can be negative to round the integral part.
    /// (This is identical to MySQL/Oracle's behavior)
    /// If you want to not update the current value,
    /// use round_to() method.
    /// NOTE: current round mode is always RoundHalfUp, which is
    /// the only behavior of MySQL.
    #[inline]
    pub fn round(&mut self, frac: isize) {
        let this_frac = self.frac();
        let intg_units = self.intg_units();
        let mut frac_units = self.frac_units();
        if frac >= this_frac as isize {
            // round precision is greater than current decimal's precision
            let round_frac_units = get_units(frac as i8);
            if round_frac_units == frac_units {
                // same number of frac units
                self.frac = frac as i8;
                return;
            }
            // copy with offset
            self.lsu
                .copy_within(..frac_units + intg_units, round_frac_units - frac_units);
            self.frac = frac as i8;
            return;
        }
        let trunc: isize;
        let round_frac: i8;
        let mut carry: i32 = 0;
        if frac < 0 {
            // round integral part
            let intg = -frac;
            if intg + 1 > self.intg() as isize {
                self.set_zero();
                return;
            }
            self.lsu.copy_within(frac_units..frac_units + intg_units, 0);
            for v in &mut self.lsu[intg_units..frac_units + intg_units] {
                *v = 0; // reset higher units to zero
            }
            frac_units = 0;
            round_frac = 0;
            trunc = intg;
        } else {
            round_frac = frac as i8;
            let mut round_frac_units = get_units(frac as i8 + 1); // how many units we need to keep for rounding
            if round_frac_units < frac_units {
                // this decimal has more frac units, drop them
                if mod9(frac as i8) == 0 {
                    // at edge of one unit
                    if self.lsu[frac_units - round_frac_units] >= HALF_UNIT as i32 {
                        carry = 1;
                    }
                    round_frac_units -= 1; // drop checked unit because carry is analyzed
                }
                self.lsu[..].copy_within(frac_units - round_frac_units..frac_units + intg_units, 0);
                for v in &mut self.lsu[round_frac_units + intg_units..frac_units + intg_units] {
                    *v = 0; // reset higher units to zero
                }
                frac_units = round_frac_units;
            } else {
                debug_assert!(round_frac_units == frac_units);
                if mod9(frac as i8) == 0 {
                    // at edge of one unit
                    if self.lsu[0] >= HALF_UNIT as i32 {
                        // check rounding on least significant unit
                        carry = 1;
                    }
                    self.lsu[..].copy_within(1..frac_units + intg_units, 0); // drop checked unit
                    self.lsu[frac_units + intg_units - 1] = 0;
                    frac_units -= 1;
                }
            }
            trunc = (frac_units * DIGITS_PER_UNIT) as isize - frac; // fractional digits to remove
        }
        round_half_up(self, intg_units, frac_units, trunc, round_frac, carry);
    }

    #[inline]
    pub fn round_to(&self, res: &mut FixedDecimal, frac: isize) {
        res.set_zero(); // always clear result first
        res.intg = self.intg; // copy intg with sign for future check
        let this_frac = self.frac();
        let intg_units = self.intg_units();
        let mut frac_units = self.frac_units();
        if frac >= this_frac as isize {
            // round precision is larger than or equal to current decimal's precision
            let round_frac_units = get_units(frac as i8);
            if round_frac_units == frac_units {
                // same number of frac units
                res.intg = self.intg;
                res.frac = frac as i8;
                res.lsu = self.lsu;
                return;
            }
            // copy with offset
            res.lsu[round_frac_units - frac_units..round_frac_units + intg_units]
                .copy_from_slice(&self.lsu[..frac_units + intg_units]);
            res.intg = self.intg;
            res.frac = frac as i8;
            return;
        }
        let trunc: isize;
        let round_frac: i8;
        let mut carry: i32 = 0;
        if frac < 0 {
            // round integral part
            let intg = -frac;
            if intg + 1 > self.intg() as isize {
                res.set_zero();
                return;
            }
            res.lsu[..intg_units].copy_from_slice(&self.lsu[frac_units..frac_units + intg_units]); // remove fractional units
            frac_units = 0;
            round_frac = 0;
            trunc = intg;
        } else {
            round_frac = frac as i8;
            let mut round_frac_units = get_units(frac as i8 + 1); // how many units we need to keep for rounding
            if round_frac_units < frac_units {
                // this decimal has more frac units, drop them
                if mod9(frac as i8) == 0 {
                    // at edge of one unit
                    if self.lsu[frac_units - round_frac_units] >= HALF_UNIT as i32 {
                        carry = 1;
                    }
                    round_frac_units -= 1; // drop checked unit because carry is analyzed
                }
                res.lsu[..round_frac_units + intg_units].copy_from_slice(
                    &self.lsu[frac_units - round_frac_units..frac_units + intg_units],
                );
                frac_units = round_frac_units;
            } else {
                debug_assert!(round_frac_units == frac_units);
                if mod9(frac as i8) == 0 {
                    // at edge of one unit
                    if self.lsu[0] >= HALF_UNIT as i32 {
                        // check rounding on least significant unit
                        carry = 1;
                    }
                    res.lsu[..frac_units + intg_units - 1]
                        .copy_from_slice(&self.lsu[1..frac_units + intg_units]); // drop checked unit
                    frac_units -= 1;
                } else {
                    res.lsu = self.lsu;
                }
            }
            trunc = (frac_units * DIGITS_PER_UNIT) as isize - frac; // fractional digits to remove
        }
        round_half_up(res, intg_units, frac_units, trunc, round_frac, carry);
    }
}

fn round_half_up(
    fd: &mut FixedDecimal,
    mut intg_units: usize,
    frac_units: usize,
    trunc: isize,
    round_frac: i8,
    mut carry: i32,
) {
    let mut round_idx = div9(trunc as i8);
    let round_pos = mod9(trunc as i8);
    let clear_idx = round_idx;
    if round_pos != 0 {
        let mut u = fd.lsu[round_idx];
        match round_pos {
            // unroll the code to optimize arithmetic with const values
            1 => {
                let r = u % 10;
                u -= r;
                if r >= 5 {
                    // round up
                    u += 10;
                }
            }
            2 => {
                let r = u % 100;
                u -= r;
                if r >= 50 {
                    // round up
                    u += 100;
                }
            }
            3 => {
                let r = u % 1000;
                u -= r;
                if r >= 500 {
                    u += 1000;
                }
            }
            4 => {
                let r = u % 10_000;
                u -= r;
                if r >= 5_000 {
                    u += 10_000;
                }
            }
            5 => {
                let r = u % 100_000;
                u -= r;
                if r >= 50_000 {
                    u += 100_000;
                }
            }
            6 => {
                let r = u % 1_000_000;
                u -= r;
                if r >= 500_000 {
                    u += 1_000_000;
                }
            }
            7 => {
                let r = u % 10_000_000;
                u -= r;
                if r >= 5_000_000 {
                    u += 10_000_000;
                }
            }
            8 => {
                let r = u % 100_000_000;
                u -= r;
                if r >= 50_000_000 {
                    u += 100_000_000;
                }
            }
            _ => unreachable!("round position unreachable"),
        }
        if u >= UNIT as i32 {
            // carry to higher unit
            u -= UNIT as i32;
            carry = 1;
        }
        fd.lsu[round_idx] = u;
        round_idx += 1;
    }
    let end_idx = intg_units + frac_units;
    if carry > 0 {
        for u in &mut fd.lsu[round_idx..end_idx] {
            let (r0, c0) = add_with_carry(*u, 0, carry);
            *u = r0;
            carry = c0;
            if carry == 0 {
                break;
            }
        }
        if carry > 0 {
            fd.lsu[end_idx] = 1;
            intg_units += 1;
        }
    }
    // clear all units below
    fd.lsu[..clear_idx].fill(0);
    let neg = fd.is_neg();
    fd.intg = (intg_units * DIGITS_PER_UNIT) as i8;
    fd.frac = round_frac;
    if neg {
        fd.set_neg();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_round() {
        let mut fd1 = FixedDecimal::zero();
        let mut fd2 = FixedDecimal::zero();
        for (input, frac, expected) in vec![
            ("0", 0, "0"),
            ("0", 1, "0.0"),
            ("1", -1, "0"),
            ("1", -2, "0"),
            ("1.00", 1, "1.0"),
            ("1.01", 1, "1.0"),
            ("1.02", 1, "1.0"),
            ("1.03", 1, "1.0"),
            ("1.04", 1, "1.0"),
            ("1.05", 1, "1.1"),
            ("1.06", 1, "1.1"),
            ("1.07", 1, "1.1"),
            ("1.08", 1, "1.1"),
            ("1.09", 1, "1.1"),
            ("1.050", 1, "1.1"),
            ("1.051", 1, "1.1"),
            ("1.052", 1, "1.1"),
            ("1.053", 1, "1.1"),
            ("1.054", 1, "1.1"),
            ("1.055", 1, "1.1"),
            ("1.056", 1, "1.1"),
            ("1.057", 1, "1.1"),
            ("1.058", 1, "1.1"),
            ("1.059", 1, "1.1"),
            ("1.040", 1, "1.0"),
            ("1.041", 1, "1.0"),
            ("1.042", 1, "1.0"),
            ("1.043", 1, "1.0"),
            ("1.044", 1, "1.0"),
            ("1.045", 1, "1.0"),
            ("1.046", 1, "1.0"),
            ("1.047", 1, "1.0"),
            ("1.048", 1, "1.0"),
            ("1.049", 1, "1.0"),
            ("1.0000000000", 9, "1.000000000"),
            ("1.0000000001", 9, "1.000000000"),
            ("1.0000000002", 9, "1.000000000"),
            ("1.0000000003", 9, "1.000000000"),
            ("1.0000000004", 9, "1.000000000"),
            ("1.0000000005", 9, "1.000000001"),
            ("1.0000000006", 9, "1.000000001"),
            ("1.0000000007", 9, "1.000000001"),
            ("1.0000000008", 9, "1.000000001"),
            ("1.0000000009", 9, "1.000000001"),
            ("1.0000000090", 9, "1.000000009"),
            ("1.0000000091", 9, "1.000000009"),
            ("1.0000000092", 9, "1.000000009"),
            ("1.0000000093", 9, "1.000000009"),
            ("1.0000000094", 9, "1.000000009"),
            ("1.0000000095", 9, "1.000000010"),
            ("1.0000000096", 9, "1.000000010"),
            ("1.0000000097", 9, "1.000000010"),
            ("1.0000000098", 9, "1.000000010"),
            ("1.0000000099", 9, "1.000000010"),
            ("999999999.0", 0, "999999999"),
            ("999999999.1", 0, "999999999"),
            ("999999999.2", 0, "999999999"),
            ("999999999.3", 0, "999999999"),
            ("999999999.4", 0, "999999999"),
            ("999999999.5", 0, "1000000000"),
            ("999999999.6", 0, "1000000000"),
            ("999999999.7", 0, "1000000000"),
            ("999999999.8", 0, "1000000000"),
            ("999999999.9", 0, "1000000000"),
            ("999999999.99990", 4, "999999999.9999"),
            ("999999999.99991", 4, "999999999.9999"),
            ("999999999.99992", 4, "999999999.9999"),
            ("999999999.99993", 4, "999999999.9999"),
            ("999999999.99994", 4, "999999999.9999"),
            ("999999999.99995", 4, "1000000000.0000"),
            ("999999999.99996", 4, "1000000000.0000"),
            ("999999999.99997", 4, "1000000000.0000"),
            ("999999999.99998", 4, "1000000000.0000"),
            ("999999999.99999", 4, "1000000000.0000"),
            ("999999999999999999.0", 0, "999999999999999999"),
            ("999999999999999999.1", 0, "999999999999999999"),
            ("999999999999999999.2", 0, "999999999999999999"),
            ("999999999999999999.3", 0, "999999999999999999"),
            ("999999999999999999.4", 0, "999999999999999999"),
            ("999999999999999999.5", 0, "1000000000000000000"),
            ("999999999999999999.6", 0, "1000000000000000000"),
            ("999999999999999999.7", 0, "1000000000000000000"),
            ("999999999999999999.8", 0, "1000000000000000000"),
            ("999999999999999999.9", 0, "1000000000000000000"),
            ("999999999999999999.90", 1, "999999999999999999.9"),
            ("999999999999999999.91", 1, "999999999999999999.9"),
            ("999999999999999999.92", 1, "999999999999999999.9"),
            ("999999999999999999.93", 1, "999999999999999999.9"),
            ("999999999999999999.94", 1, "999999999999999999.9"),
            ("999999999999999999.95", 1, "1000000000000000000.0"),
            ("999999999999999999.96", 1, "1000000000000000000.0"),
            ("999999999999999999.97", 1, "1000000000000000000.0"),
            ("999999999999999999.98", 1, "1000000000000000000.0"),
            ("999999999999999999.99", 1, "1000000000000000000.0"),
            ("0.9876543210", 10, "0.9876543210"),
            ("0.9876543210", 9, "0.987654321"),
            ("0.9876543210", 8, "0.98765432"),
            ("0.9876543210", 7, "0.9876543"),
            ("0.9876543210", 6, "0.987654"),
            ("0.9876543210", 5, "0.98765"),
            ("0.9876543210", 4, "0.9877"),
            ("0.9876543210", 3, "0.988"),
            ("0.9876543210", 2, "0.99"),
            ("0.9876543210", 1, "1.0"),
            ("0.9876543210", 0, "1"),
            ("-1.9876543210", 0, "-2"),
            ("-1.9876543210", 1, "-2.0"),
            ("123456789123456789", -1, "123456789123456790"),
            ("123456789123456789", -2, "123456789123456800"),
            ("123456789123456789", -3, "123456789123457000"),
            ("123456789123456789", -4, "123456789123460000"),
            ("123456789123456789", -5, "123456789123500000"),
            ("123456789123456789", -6, "123456789123000000"),
            ("123456789123456789", -7, "123456789120000000"),
            ("123456789123456789", -8, "123456789100000000"),
            ("123456789123456789", -9, "123456789000000000"),
            ("999999999999999999", -1, "1000000000000000000"),
            ("0.999999999", 9, "0.999999999"),
            ("0.999999999", 8, "1.00000000"),
            ("0.999999999", 7, "1.0000000"),
            ("0.999999999", 6, "1.000000"),
            ("0.999999999", 5, "1.00000"),
            ("0.999999999", 4, "1.0000"),
            ("0.999999999", 3, "1.000"),
            ("0.999999999", 2, "1.00"),
            ("0.999999999", 1, "1.0"),
            ("0.999999999", 0, "1"),
        ] {
            assert!(fd1.from_ascii_str(input, true).is_ok());
            fd1.round_to(&mut fd2, frac);
            let actual = fd2.to_string(-1);
            println!("({}).round({}) = {}", fd1.to_string(-1), frac, actual);
            assert_eq!(&actual, expected);
            fd1.round(frac);
            let actual = fd1.to_string(-1);
            println!("self.round({}) = {}", frac, fd1.to_string(-1));
            assert_eq!(&actual, expected);
        }
    }
}
