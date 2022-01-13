use crate::decimal::FixedDecimal;
use std::cmp::Ordering;

impl PartialEq for FixedDecimal {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for FixedDecimal {}

impl PartialOrd for FixedDecimal {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FixedDecimal {
    fn cmp(&self, other: &Self) -> Ordering {
        let lneg = self.is_neg();
        let rneg = other.is_neg();
        if lneg {
            if rneg {
                return cmp_abs(other, self);
            }
            return Ordering::Less;
        }
        if rneg {
            return Ordering::Greater;
        }
        cmp_abs(self, other)
    }
}

#[inline]
pub(crate) fn cmp_abs(lhs: &FixedDecimal, rhs: &FixedDecimal) -> Ordering {
    cmp_abs_lsu(
        lhs.intg_units() as isize, lhs.frac_units() as isize, &lhs.lsu[..], 
        rhs.intg_units() as isize, rhs.frac_units() as isize, &rhs.lsu[..])
}

pub(crate) fn cmp_abs_lsu(mut liu: isize, lfu: isize, llsu: &[i32], mut riu: isize, rfu: isize, rlsu: &[i32]) -> Ordering {
    let mut i = liu + lfu - 1;
    let mut j = riu + rfu - 1;
    while liu > 0 && liu > riu { // lhs has more integral units
        if llsu[i as usize] > 0 {
            return Ordering::Greater;
        }
        liu -= 1;
        i -= 1;
    }
    while riu > 0 && riu > liu { // rhs has more integral units
        if rlsu[j as usize] > 0 {
            return Ordering::Less;
        }
        riu -= 1;
        j -= 1;
    }
    // both have identical number of integral units
    while i >= 0 && j >= 0 {
        let lv = llsu[i as usize];
        let rv = rlsu[j as usize];
        if lv > rv {
            return Ordering::Greater;
        }
        if lv < rv {
            return Ordering::Less;
        }
        i -= 1;
        j -= 1;
    }
    while i >= 0 { // lhs still has units
        if llsu[i as usize] > 0 {
            return Ordering::Greater;
        }
        i -= 1;
    }
    while j >= 0 {
        if rlsu[j as usize] > 0 {
            return Ordering::Less;
        }
        j-= 1;
    }
    // no units to compare so they are equal
    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cmp() {
        let mut fd1 = FixedDecimal::zero();
        let mut fd2 = FixedDecimal::zero();
        for (input1, input2, expected) in vec![
            ("0", "0", Ordering::Equal),
            ("0", "1", Ordering::Less),
            ("1", "0", Ordering::Greater),
            ("-1", "-1", Ordering::Equal),
            ("1", "-1", Ordering::Greater),
            ("-1", "1", Ordering::Less),
            ("-1", "-2", Ordering::Greater),
            ("1", "1", Ordering::Equal),
            ("2", "1", Ordering::Greater),
            ("1", "2", Ordering::Less),
            ("1.0", "1", Ordering::Equal),
            ("1", "1.0", Ordering::Equal),
            ("1.000", "1.00", Ordering::Equal),
            ("1.000000000000", "1.00000", Ordering::Equal),
            ("1.01", "1", Ordering::Greater),
            ("1", "1.01", Ordering::Less),
            ("1.02", "1.01", Ordering::Greater),
            ("1.01", "1.02", Ordering::Less),
            ("1000000000", "999999999", Ordering::Greater),
            ("999999999", "1000000000", Ordering::Less),
            ("1.0000000000000000000000000001", "1", Ordering::Greater),
            ("1", "1.0000000000000000000000000001", Ordering::Less),
            ("1.0000000010000000000000000001", "1.0000000010000000000000000001", Ordering::Equal),
            ("1.0000000010000000000000000001", "1.0000000010000000000000000002", Ordering::Less),
            ("1.0000000010000000000000000002", "1.0000000010000000000000000001", Ordering::Greater),
            ("1.0000000010000000010000000001", "1.0000000010000000010000000001", Ordering::Equal),
            ("1.0000000010000000010000000001", "1.0000000010000000010000000002", Ordering::Less),
            ("1.0000000010000000010000000002", "1.0000000010000000010000000001", Ordering::Greater),
            ("100000000000000000000000000", "100000000000000000000000000", Ordering::Equal),
            ("100000000000000000000000001", "100000000000000000000000000", Ordering::Greater),
            ("100000000000000000000000000", "100000000000000000000000001", Ordering::Less),
        ] {
            assert!(fd1.from_ascii_str(input1, true).is_ok());
            assert!(fd2.from_ascii_str(input2, true).is_ok());
            let actual = fd1.cmp(&fd2);
            assert_eq!(expected, actual);
        }
    }
}