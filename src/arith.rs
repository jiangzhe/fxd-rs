use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Rem, Sub};

use crate::cmp::cmp_abs_lsu;
use crate::decimal::{FixedDecimal, DIGITS_PER_UNIT, MAX_FRAC_UNITS, MAX_UNITS, UNIT};
use crate::error::{Error, Result};
use crate::utils::{get_units, unit_leading_zeroes, units_greater_equal};

const DIV_INCR_FRAC: usize = 4;

impl FixedDecimal {
    #[inline]
    pub fn add_to(lhs: &Self, rhs: &Self, res: &mut Self) -> Result<()> {
        if lhs.is_zero() {
            *res = rhs.clone();
            return Ok(());
        }
        if rhs.is_zero() {
            *res = lhs.clone();
            return Ok(());
        }
        let lneg = lhs.is_neg();
        let rneg = rhs.is_neg();
        if lneg == rneg {
            add_abs(lhs, rhs, res)?;
            if lneg {
                res.set_neg_and_check_zero();
            }
            return Ok(());
        }
        let sub_neg = sub_abs(lhs, rhs, res)?;
        if sub_neg != lneg {
            res.set_neg_and_check_zero();
        }
        Ok(())
    }

    #[inline]
    pub fn sub_to(lhs: &Self, rhs: &Self, res: &mut Self) -> Result<()> {
        if lhs.is_zero() {
            *res = rhs.clone();
            if res.is_neg() {
                res.set_pos();
            } else {
                res.set_neg_and_check_zero();
            }
            return Ok(());
        }
        if rhs.is_zero() {
            *res = lhs.clone();
            return Ok(());
        }
        let lneg = lhs.is_neg();
        let rneg = rhs.is_neg();
        if lneg != rneg {
            add_abs(lhs, rhs, res)?;
            if lneg {
                res.set_neg_and_check_zero();
            }
            return Ok(());
        }
        let sub_neg = sub_abs(lhs, rhs, res)?;
        if sub_neg != lneg {
            res.set_neg_and_check_zero();
        }
        Ok(())
    }

    #[inline]
    pub fn mul_to(lhs: &Self, rhs: &Self, res: &mut Self) -> Result<()> {
        if lhs.is_zero() || rhs.is_zero() {
            res.set_zero();
            return Ok(());
        }
        let res_neg = lhs.is_neg() != rhs.is_neg();
        mul_abs(lhs, rhs, res)?;
        if res_neg {
            res.set_neg_and_check_zero();
        }
        Ok(())
    }

    #[inline]
    pub fn div_to(lhs: &Self, rhs: &Self, res: &mut Self, incr_frac: usize) -> Result<()> {
        let res_neg = lhs.is_neg() != rhs.is_neg();
        div_abs(lhs, rhs, res, incr_frac as isize)?;
        if res_neg {
            res.set_neg_and_check_zero();
        }
        Ok(())
    }

    #[inline]
    pub fn rem_to(lhs: &Self, rhs: &Self, res: &mut Self) -> Result<()> {
        let res_neg = lhs.is_neg();
        rem_abs(lhs, rhs, res)?;
        if res_neg {
            res.set_neg_and_check_zero();
        }
        Ok(())
    }
}

impl Add for FixedDecimal {
    type Output = FixedDecimal;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        let mut res = FixedDecimal {
            intg: 0,
            frac: 0,
            lsu: [0; MAX_UNITS],
        };
        if let Err(e) = Self::add_to(&self, &rhs, &mut res) {
            panic!("addition failed: {}", e);
        }
        res
    }
}

impl Sub for FixedDecimal {
    type Output = FixedDecimal;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        let mut res = FixedDecimal {
            intg: 0,
            frac: 0,
            lsu: [0; MAX_UNITS],
        };
        if let Err(e) = Self::sub_to(&self, &rhs, &mut res) {
            panic!("subtraction failed: {}", e);
        }
        res
    }
}

impl Mul for FixedDecimal {
    type Output = FixedDecimal;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = FixedDecimal {
            intg: 0,
            frac: 0,
            lsu: [0; MAX_UNITS],
        };
        if let Err(e) = Self::mul_to(&self, &rhs, &mut res) {
            panic!("multiplication failed: {}", e);
        }
        res
    }
}

impl Div for FixedDecimal {
    type Output = FixedDecimal;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        let mut res = FixedDecimal {
            intg: 0,
            frac: 0,
            lsu: [0; MAX_UNITS],
        };
        if let Err(e) = Self::div_to(&self, &rhs, &mut res, DIV_INCR_FRAC) {
            panic!("division failed: {}", e);
        }
        res
    }
}

impl Rem for FixedDecimal {
    type Output = FixedDecimal;

    #[inline]
    fn rem(self, rhs: Self) -> Self::Output {
        let mut res = FixedDecimal {
            intg: 0,
            frac: 0,
            lsu: [0; MAX_UNITS],
        };
        if let Err(e) = Self::rem_to(&self, &rhs, &mut res) {
            panic!("modulo failed: {}", e);
        }
        res
    }
}

/// sums two decimals' absolute values.
/// Separate units into 3 segments, intg, common, frac
/// | lhs  |  xxxx  |  xxxx.xxxx  |        |
/// | rhs  |        |  yyyy.yyyy  |  yyyy  |
/// |------|--------|-------------|--------|
/// |      |  intg  |   common    |  frac  |
///
/// frac can be directly added to result.
/// common uses normal addition, taking care of carry.
/// intg is added with 0 and carry.
fn add_abs(lhs: &FixedDecimal, rhs: &FixedDecimal, res: &mut FixedDecimal) -> Result<()> {
    res.set_zero(); // always clear result first
    let liu = lhs.intg_units();
    let lfu = lhs.frac_units();
    let riu = rhs.intg_units();
    let rfu = rhs.frac_units();
    let mut frac_unit_diff = lfu as isize - rfu as isize;
    let mut lhs_idx: isize = 0;
    let mut rhs_idx: isize = 0;
    let mut res_idx: isize = 0;
    match frac_unit_diff.cmp(&0) {
        Ordering::Greater => {
            res.lsu[..frac_unit_diff as usize].copy_from_slice(&lhs.lsu[..frac_unit_diff as usize]);
            lhs_idx = frac_unit_diff;
            res_idx = frac_unit_diff;
        }
        Ordering::Less => {
            frac_unit_diff = -frac_unit_diff;
            res.lsu[..frac_unit_diff as usize].copy_from_slice(&rhs.lsu[..frac_unit_diff as usize]);
            rhs_idx = frac_unit_diff;
            res_idx = frac_unit_diff;
        }
        Ordering::Equal => (),
    }
    let mut stop = res_idx + liu.min(riu) as isize + lfu.min(rfu) as isize;
    let mut carry: i32 = 0;
    while res_idx < stop {
        let (r0, c0) = add_with_carry(lhs.lsu[lhs_idx as usize], rhs.lsu[rhs_idx as usize], carry);
        res.lsu[res_idx as usize] = r0;
        carry = c0;
        lhs_idx += 1;
        rhs_idx += 1;
        res_idx += 1;
    }
    let intg_unit_diff = liu as isize - riu as isize;
    if intg_unit_diff > 0 {
        // lhs has more intg
        stop = lhs_idx + intg_unit_diff;
        while lhs_idx < stop {
            let (r0, c0) = add_with_carry(lhs.lsu[lhs_idx as usize], 0, carry);
            res.lsu[res_idx as usize] = r0;
            carry = c0;
            res_idx += 1;
            lhs_idx += 1;
        }
        if carry > 0 {
            res.lsu[res_idx as usize] = carry;
        }
    } else if intg_unit_diff < 0 {
        // rhs has more intg
        stop = rhs_idx - intg_unit_diff;
        while rhs_idx < stop {
            let (r0, c0) = add_with_carry(rhs.lsu[rhs_idx as usize], 0, carry);
            res.lsu[res_idx as usize] = r0;
            carry = c0;
            res_idx += 1;
            rhs_idx += 1;
        }
        if carry > 0 {
            res.lsu[res_idx as usize] = carry;
        }
    } else if carry != 0 {
        // no more intg but carry is non-zero
        res.lsu[res_idx as usize] = carry;
    }

    let mut res_intg_non_zero = liu.max(riu) as isize + carry as isize;
    res.frac = lhs.frac.max(rhs.frac);
    let res_frac_units = res.frac_units() as isize;
    while res_intg_non_zero >= 0 {
        if res.lsu[(res_intg_non_zero + res_frac_units) as usize] > 0 {
            break;
        }
        res_intg_non_zero -= 1;
    }
    if res_intg_non_zero >= 0 {
        res.intg = ((res_intg_non_zero + 1) * DIGITS_PER_UNIT as isize) as i8;
    } else {
        res.intg = 0;
    }
    Ok(())
}

/// differs two decimals' absolute values, returns the negative flag.
/// Similar to addAbs(), but if lhs is smaller than rhs, we will have
/// borrow=-1 at end. Then we can traverse all units and apply subtraction
/// again.
///
/// Separate units into 3 segments, intgSeg, commonSeg, fracSeg
/// | lhs  |  xxxx  |  xxxx.xxxx  |        |
/// | rhs  |        |  yyyy.yyyy  |  yyyy  |
/// |------|--------|-------------|--------|
/// |      |  intg  |   common    |  frac  |
fn sub_abs(lhs: &FixedDecimal, rhs: &FixedDecimal, res: &mut FixedDecimal) -> Result<bool> {
    res.set_zero(); // always clear result first
    let liu = lhs.intg_units();
    let lfu = lhs.frac_units();
    let riu = rhs.intg_units();
    let rfu = rhs.frac_units();
    let mut frac_unit_diff = lfu as isize - rfu as isize;
    let mut lhs_idx: isize = 0;
    let mut rhs_idx: isize = 0;
    let mut res_idx: isize = 0;
    let mut borrow: i32 = 0;
    match frac_unit_diff.cmp(&0) {
        // lhs has more frac
        Ordering::Greater => {
            res.lsu[..frac_unit_diff as usize].copy_from_slice(&lhs.lsu[..frac_unit_diff as usize]);
            lhs_idx = frac_unit_diff;
            res_idx = frac_unit_diff;
        }
        Ordering::Less => {
            frac_unit_diff = -frac_unit_diff;
            for (rst, rv) in res
                .lsu
                .iter_mut()
                .zip(rhs.lsu.iter())
                .take(frac_unit_diff as usize)
            {
                let (r0, b0) = sub_with_borrow(0, *rv, borrow);
                *rst = r0;
                borrow = b0;
            }
            rhs_idx = frac_unit_diff;
            res_idx = frac_unit_diff;
        }
        Ordering::Equal => (),
    }
    let mut stop = res_idx + liu.min(riu) as isize + lfu.min(rfu) as isize;
    while res_idx < stop {
        let (r0, b0) =
            sub_with_borrow(lhs.lsu[lhs_idx as usize], rhs.lsu[rhs_idx as usize], borrow);
        res.lsu[res_idx as usize] = r0;
        borrow = b0;
        res_idx += 1;
        lhs_idx += 1;
        rhs_idx += 1;
    }
    let intg_unit_diff = liu as isize - riu as isize;
    match intg_unit_diff.cmp(&0) {
        Ordering::Greater => {
            stop = lhs_idx + intg_unit_diff;
            while lhs_idx < stop {
                let (r0, b0) = sub_with_borrow(lhs.lsu[lhs_idx as usize], 0, borrow);
                res.lsu[res_idx as usize] = r0;
                borrow = b0;
                res_idx += 1;
                lhs_idx += 1;
            }
        }
        Ordering::Less => {
            stop = rhs_idx - intg_unit_diff;
            while rhs_idx < stop {
                let (r0, b0) = sub_with_borrow(0, rhs.lsu[rhs_idx as usize], borrow);
                res.lsu[res_idx as usize] = r0;
                borrow = b0;
                res_idx += 1;
                rhs_idx += 1;
            }
        }
        Ordering::Equal => (),
    }
    let neg = borrow == -1;
    if neg {
        borrow = 0;
        for i in 0..res_idx {
            let (r0, b0) = sub_with_borrow(0, res.lsu[i as usize], borrow);
            res.lsu[i as usize] = r0;
            borrow = b0;
        }
    }
    let mut res_intg_non_zero = liu.max(riu) as isize;
    res.frac = lhs.frac.max(rhs.frac);
    let result_frac_units = res.frac_units();
    while res_intg_non_zero >= 0 {
        if res.lsu[(res_intg_non_zero + result_frac_units as isize) as usize] > 0 {
            break;
        }
        res_intg_non_zero -= 1;
    }
    if res_intg_non_zero >= 0 {
        res.intg = ((res_intg_non_zero + 1) * DIGITS_PER_UNIT as isize) as i8;
    } else {
        res.intg = 0;
    }
    Ok(neg)
}

/// multiplies two decimals' absolute values.
/// The result frac precision is extended to sum of both precisions. until
/// reaching the limitation of MAX_UNITS or MAX_FRAC_UNITS.
fn mul_abs(lhs: &FixedDecimal, rhs: &FixedDecimal, res: &mut FixedDecimal) -> Result<()> {
    res.set_zero(); // always clear result first
    let res_intg_digits = lhs.intg() + rhs.intg();
    let res_intg_units = get_units(res_intg_digits);
    let res_frac_digits = lhs.frac() + rhs.frac();
    let mut res_frac_units = get_units(res_frac_digits);
    if res_intg_units > MAX_UNITS {
        // integral overflow
        return Err(Error::Overflow);
    }
    if res_intg_units + res_frac_units > MAX_UNITS {
        // integral+fractional overflow
        res_frac_units = MAX_UNITS - res_intg_units; // fractional truncation required
    }
    if res_frac_units > MAX_FRAC_UNITS {
        // exceeds maximum fractional units
        res_frac_units = MAX_FRAC_UNITS;
    }
    let liu = lhs.intg_units();
    let lfu = lhs.frac_units();
    let riu = rhs.intg_units();
    let rfu = rhs.frac_units();
    // because result fractional part may be truncated, we need to calculate
    // how many units has to be shifted and can be ignored in calculation.
    // 1. If lhsIdx + rhsIdx - shiftUnits < -1, we can ignore the result
    //    of left.lsu[lhsIdx]*right.lsu[rhsIdx].
    // 2. If lhsIdx + rhsIdx - shiftUnits = -1, we need to take care of the
    //    carry and add to result.lsu[0].
    // 3. If lhsIdx + rhsIdx - shiftUnits >= 0, follow normal calcuation.
    let shift_units = (lfu + rfu - res_frac_units) as isize;
    let mut carry: i64 = 0;
    // let mut res_idx: isize = 0;
    for (rhs_idx, rv) in rhs.lsu[..(riu + rfu) as usize].iter().enumerate() {
        let mut res_idx = -1;
        for (lhs_idx, lv) in lhs.lsu[..(liu + lfu) as usize].iter().enumerate() {
            res_idx = (lhs_idx + rhs_idx) as isize - shift_units;
            if res_idx < -1 {
                continue;
            }
            if res_idx == -1 {
                let v = (*lv as i64) * (*rv as i64);
                if v < UNIT as i64 {
                    continue;
                }
                carry = v / UNIT as i64; // we only need the carry for lsu of result
                continue;
            }
            // calculate product and sum with previous result and carry
            let v = *lv as i64 * *rv as i64 + res.lsu[res_idx as usize] as i64 + carry;
            carry = v / UNIT as i64;
            res.lsu[res_idx as usize] = (v - carry * UNIT as i64) as i32;
        }
        debug_assert!(res_idx >= 0);
        if res_idx + 1 < MAX_UNITS as isize {
            res.lsu[(res_idx + 1) as usize] = carry as i32;
        } else if carry > 0 {
            return Err(Error::Overflow);
        }
        carry = 0
    }
    res.frac = res_frac_digits.min((res_frac_units * DIGITS_PER_UNIT) as i8);
    res.intg = (res_intg_units * DIGITS_PER_UNIT) as i8;
    Ok(())
}

/// Divides two decimals' absolute values.
/// It's implementation of Knuth's Algorithm 4.3.1 D, with support on fractional numbers.
#[allow(clippy::many_single_char_names)]
fn div_abs(
    lhs: &FixedDecimal,
    rhs: &FixedDecimal,
    res: &mut FixedDecimal,
    mut incr_frac: isize,
) -> Result<()> {
    res.set_zero(); // always clear result first
    let lhs_intg = lhs.intg();
    let liu = get_units(lhs_intg);
    let lhs_frac = lhs.frac();
    let lfu = get_units(lhs_frac);
    let lhs_ext_frac = (lfu * DIGITS_PER_UNIT) as isize; // extended frac with unit size
    let rhs_intg = rhs.intg();
    let riu = get_units(rhs_intg);
    let rhs_frac = rhs.frac();
    let rfu = get_units(rhs_frac);
    let rhs_ext_frac = (rfu * DIGITS_PER_UNIT) as isize; // extended frac with unit size
                                                         // leading non-zero unit of lhs and rhs
    let mut rhs_non_zero = (riu + rfu) as isize - 1;
    while rhs_non_zero >= 0 {
        if rhs.lsu[rhs_non_zero as usize] > 0 {
            break;
        }
        rhs_non_zero -= 1;
    }
    if rhs_non_zero < 0 {
        // divider is zero
        return Err(Error::DivisionByZero);
    }

    // check and remove leading zeroes in lhs
    let mut lhs_non_zero = (liu + lfu) as isize - 1;
    while lhs_non_zero >= 0 {
        if lhs.lsu[lhs_non_zero as usize] > 0 {
            break;
        }
        lhs_non_zero -= 1;
    }
    if lhs_non_zero < 0 {
        // dividend is zero
        return Ok(());
    }
    // digits of rhs from leading non-zero position
    let rhs_prec = (rhs_non_zero + 1) * DIGITS_PER_UNIT as isize
        - unit_leading_zeroes(rhs.lsu[rhs_non_zero as usize]) as isize;
    // digits of lhs from leading non-zero position
    let lhs_prec = (lhs_non_zero + 1) * DIGITS_PER_UNIT as isize
        - unit_leading_zeroes(lhs.lsu[lhs_non_zero as usize]) as isize;
    // because we store fractional part in units, we always extend frac precision
    // to multiple of 9. Here check if incrPrec is already covered by the auto-extension.
    // if so, reset incrPrec to 0.
    incr_frac -= lhs_ext_frac - lhs_frac as isize + rhs_ext_frac - rhs_frac as isize;
    if incr_frac < 0 {
        incr_frac = 0;
    }
    // calculate result frac units
    let mut res_frac_units = get_units((lhs_ext_frac + rhs_ext_frac + incr_frac) as i8);
    if res_frac_units > MAX_FRAC_UNITS {
        // exceeds maximum fractional digits
        res_frac_units = MAX_FRAC_UNITS;
    }
    // calculate result intg units
    let mut res_intg = (lhs_prec - lhs_ext_frac) - (rhs_prec - rhs_ext_frac);
    let mut dividend_shift: isize = 0;
    if units_greater_equal(
        &lhs.lsu[..(lhs_non_zero + 1) as usize],
        &rhs.lsu[..(rhs_non_zero + 1) as usize],
    ) {
        res_intg += 1;
    } else {
        dividend_shift = -1;
    }
    // now adjust result units based on limitation of maximum precision and determine
    // the start position of result unit.
    let res_intg_units: usize;
    let res_start_idx: isize;
    if res_intg > 0 {
        res_intg_units = get_units(res_intg as i8);
        if res_intg_units > MAX_UNITS {
            return Err(Error::Overflow);
        }
        if res_intg_units + res_frac_units > MAX_UNITS {
            res_frac_units = MAX_UNITS - res_intg_units;
        }
        res_start_idx = (res_frac_units + res_intg_units) as isize - 1;
    } else {
        res_intg_units = 0;
        let res_start_offset = get_units((1 - res_intg) as i8) as isize;
        res_start_idx = res_frac_units as isize - res_start_offset;
        res_intg = 0;
    }
    let res_units = res_intg_units + res_frac_units;

    // here we identify short/long division
    // short division means the divider only has single unit.
    if rhs_non_zero == 0 {
        // short division
        let d = rhs.lsu[0] as i64;
        let mut rem: i64 = 0;
        if dividend_shift < 0 {
            rem = lhs.lsu[lhs_non_zero as usize] as i64;
        }
        let mut i = lhs_non_zero + dividend_shift; // index of lhs
        let mut j = res_start_idx; // index of result
        while j >= 0 {
            let u = if i >= 0 {
                rem * UNIT as i64 + lhs.lsu[i as usize] as i64
            } else {
                rem * UNIT as i64
            };
            let q = u / d; // div
            rem = u - q * d; // update remainder
            res.lsu[j as usize] = q as i32; // update result
            i -= 1;
            j -= 1;
        }
        res.intg = res_intg as i8;
        res.frac = (res_frac_units * DIGITS_PER_UNIT) as i8;
        return Ok(());
    }

    // long division using Knuth's algorithm
    // D1. normalization
    let norm_factor = UNIT as i64 / (rhs.lsu[rhs_non_zero as usize] as i64 + 1);
    // normalize buf1 and buf2
    let mut buf1 = [0i32; MAX_UNITS * 2]; // store normalized lhs
    let mut buf2 = [0i32; MAX_UNITS]; // store normalized rhs
    if norm_factor == 1 {
        buf1[res_units..res_units + lhs_non_zero as usize + 1]
            .copy_from_slice(&lhs.lsu[..lhs_non_zero as usize + 1]);
        buf2[..rhs_non_zero as usize + 1].copy_from_slice(&rhs.lsu[..rhs_non_zero as usize + 1]);
    } else {
        let mut carry: i64 = 0;
        // normalize lhs into buf1
        let buf1_end = res_units + lhs_non_zero as usize + 1;
        for (b1, lv) in buf1[res_units..buf1_end]
            .iter_mut()
            .zip(lhs.lsu[..lhs_non_zero as usize + 1].iter())
        {
            let u = *lv as i64 * norm_factor + carry;
            carry = u / UNIT as i64;
            *b1 = (u - carry * UNIT as i64) as i32;
        }
        if carry > 0 {
            buf1[buf1_end] = carry as i32;
            carry = 0;
        }
        // normalize rhs into buf2
        for (b2, rv) in buf2[..rhs_non_zero as usize + 1]
            .iter_mut()
            .zip(rhs.lsu[..rhs_non_zero as usize + 1].iter())
        {
            let u = *rv as i64 * norm_factor + carry;
            carry = u / UNIT as i64;
            *b2 = (u - carry * UNIT as i64) as i32;
        }
        debug_assert!(carry == 0, "carry must be zero in divider normalization");
    }

    let d0 = buf2[rhs_non_zero as usize] as i64; // rhs most significant unit
    let d1 = buf2[rhs_non_zero as usize - 1] as i64; // rhs second significant unit
    let mut i = res_units as isize + lhs_non_zero + dividend_shift; // index of normalized lhs
    let mut j = res_start_idx; // index of result
    while j >= 0 {
        // D3. make the guess on u1
        let u0 = buf1[i as usize + 1] as i64;
        let mut u1: i64 = 0;
        if i >= 0 {
            u1 = buf1[i as usize] as i64;
        }
        let u = u0 * UNIT as i64 + u1;
        let mut qhat = u / d0;
        let mut rhat = u - qhat * d0;
        // qhat cannot be greater or equal to UNIT
        debug_assert!(qhat < UNIT as i64, "qhat must be less than UNIT");
        let mut u2: i64 = 0;
        if i > 0 {
            u2 = buf1[i as usize - 1] as i64;
        }
        while qhat * d1 > rhat * UNIT as i64 + u2 {
            // check if qhat can satisfy next unit
            qhat -= 1; // decrease qhat
            rhat += d0; // increase rhat
        }
        // D4. multiply and subtract
        let mut ms_idx: isize = i - rhs_non_zero;
        let mut k: isize = 0;
        let mut carry: i64 = 0;
        let mut borrow: i32 = 0;
        while k <= rhs_non_zero {
            let mul_v = qhat * buf2[k as usize] as i64 + carry; // mul
            carry = mul_v / UNIT as i64; // update carry
            let mul_v0 = mul_v - carry * UNIT as i64; // in current unit
            if ms_idx < 0 {
                let (_, b0) = sub_with_borrow(0, mul_v0 as i32, borrow);
                borrow = b0;
            } else {
                let (sub_v, b0) = sub_with_borrow(buf1[ms_idx as usize], mul_v0 as i32, borrow);
                buf1[ms_idx as usize] = sub_v;
                borrow = b0;
            }
            k += 1;
            ms_idx += 1;
        }
        borrow += buf1[ms_idx as usize] - carry as i32;
        if borrow == -1 {
            // qhat is larger, cannot satisfy the whole decimal
            // D6. add back (reverse subtract)
            qhat -= 1; // decrease qhat
            borrow = 0;
            for i in 0..rhs_non_zero {
                // reverse subtract
                let (r0, b0) = sub_with_borrow(0, buf1[(i + j) as usize], borrow);
                buf1[(i + j) as usize] = r0;
                borrow = b0;
            }
        } else {
            buf1[ms_idx as usize] = 0; // clear buf1 because multiply w/ subtract succeeds
        }
        res.lsu[j as usize] = qhat as i32;
        i -= 1;
        j -= 1;
    }
    res.intg = res_intg as i8;
    res.frac = (res_frac_units * DIGITS_PER_UNIT) as i8;
    Ok(())
}

/// Modulos two decimals' absolute values.
fn rem_abs(lhs: &FixedDecimal, rhs: &FixedDecimal, res: &mut FixedDecimal) -> Result<()> {
    res.set_zero(); // always clear result first
    let lhs_intg = lhs.intg();
    let liu = get_units(lhs_intg);
    let lhs_frac = lhs.frac();
    let lfu = get_units(lhs_frac);
    let rhs_intg = rhs.intg();
    let riu = get_units(rhs_intg);
    let rhs_frac = rhs.frac();
    let rfu = get_units(rhs_frac);
    // leading non-zero unit of lhs and rhs
    let mut rhs_non_zero = (riu + rfu) as isize - 1;
    while rhs_non_zero >= 0 {
        if rhs.lsu[rhs_non_zero as usize] > 0 {
            break;
        }
        rhs_non_zero -= 1;
    }
    if rhs_non_zero < 0 {
        // divider is zero
        return Err(Error::DivisionByZero);
    }
    // check and remove leading zeroes in lhs
    let mut lhs_non_zero = (liu + lfu) as isize - 1;
    while lhs_non_zero >= 0 {
        if lhs.lsu[lhs_non_zero as usize] > 0 {
            break;
        }
        lhs_non_zero -= 1;
    }
    if lhs_non_zero < 0 {
        // dividend is zero
        return Ok(());
    }

    match cmp_abs_lsu(
        liu as isize,
        lfu as isize,
        &lhs.lsu[..],
        riu as isize,
        rfu as isize,
        &rhs.lsu[..],
    ) {
        Ordering::Less => {
            if rfu > lfu {
                // rhs has higher fractional precision
                res.lsu[rfu - lfu..liu + rfu].copy_from_slice(&lhs.lsu[..lfu + liu]);
                res.frac = rhs_frac;
                res.intg = lhs_intg;
                return Ok(());
            }
            // lhs has higher or same fractional precision
            res.lsu[..lfu + liu].copy_from_slice(&lhs.lsu[..lfu + liu]);
            res.frac = lhs_frac.max(rhs_frac);
            res.intg = lhs_intg;
            return Ok(());
        }
        Ordering::Equal => return Ok(()), // lhs equals to rhs, result is zero
        Ordering::Greater => (),
    }
    // calculate remainder frac units
    let rem_frac = lhs_frac.max(rhs_frac);
    let rem_frac_units = get_units(rem_frac);
    let dividend_shift: isize = if units_greater_equal(
        &lhs.lsu[..(lhs_non_zero + 1) as usize],
        &rhs.lsu[..(rhs_non_zero + 1) as usize],
    ) {
        0
    } else {
        -1
    };

    // align frac units between lhs and rhs
    let (lhs_shift_units, rhs_shift_units) = if lfu < rfu {
        (rfu as isize - lfu as isize, 0)
    } else {
        (0, lfu as isize - rfu as isize)
    };

    if rhs_non_zero == 0 {
        // short division
        let d = rhs.lsu[0] as i64; // single divider
        let mut buf = [0i32; MAX_UNITS * 2];
        let buf_len = lhs_shift_units + lhs_non_zero + 1;
        buf[lhs_shift_units as usize..buf_len as usize]
            .copy_from_slice(&lhs.lsu[..(lhs_non_zero + 1) as usize]);
        let mut rem = 0;
        if dividend_shift < 0 {
            rem = lhs.lsu[lhs_non_zero as usize] as i64;
        }
        let stop = rhs_non_zero + rhs_shift_units;
        let mut i: isize = buf_len - 1 + dividend_shift; // i is lhs index
        while i >= stop {
            let u = rem * UNIT as i64 + buf[i as usize] as i64;
            let q = u / d; // div
            rem = u - q * d; // update remainder
            i -= 1;
        }
        let mut res_non_zero = -1;
        if rem > 0 {
            res.lsu[(i + 1) as usize] = rem as i32;
            res_non_zero = i + 1;
        }
        while i >= 0 {
            // copy rest of lhs into result
            let v = buf[i as usize];
            res.lsu[i as usize] = v;
            if v > 0 && res_non_zero < 0 {
                res_non_zero = i;
            }
            i -= 1;
        }
        if res_non_zero >= rem_frac_units as isize {
            let res_intg_units = (res_non_zero + 1 - rem_frac_units as isize) as i8;
            res.intg = (res_intg_units * DIGITS_PER_UNIT as i8)
                - unit_leading_zeroes(res.lsu[(res_intg_units + rem_frac_units as i8 - 1) as usize])
                    as i8;
        } else {
            res.intg = 0;
        }
        res.frac = rem_frac;
        return Ok(());
    }
    let buf1_len = lhs_non_zero + 1 + lhs_shift_units;
    let buf2_len = rhs_non_zero + 1 + rhs_shift_units;
    // D1. normalization
    let norm_factor = UNIT as i64 / (rhs.lsu[rhs_non_zero as usize] as i64 + 1);
    // normalize buf1 and buf2
    let mut buf1 = [0i32; MAX_UNITS * 2];
    let mut buf2 = [0i32; MAX_UNITS * 2];
    if norm_factor == 1 {
        buf1[lhs_shift_units as usize..buf1_len as usize]
            .copy_from_slice(&lhs.lsu[..(lhs_non_zero + 1) as usize]);
        buf2[rhs_shift_units as usize..buf2_len as usize]
            .copy_from_slice(&rhs.lsu[..(rhs_non_zero + 1) as usize]);
    } else {
        let mut carry: i64 = 0;
        let mut i = 0;
        while i <= lhs_non_zero {
            let u = lhs.lsu[i as usize] as i64 * norm_factor + carry;
            carry = u / UNIT as i64; // update carry
            buf1[(i + lhs_shift_units) as usize] = (u - carry * UNIT as i64) as i32;
            i += 1;
        }
        if carry > 0 {
            buf1[buf1_len as usize] = carry as i32;
            carry = 0;
        }
        // normalize rhs into buf2
        i = 0;
        while i <= rhs_non_zero {
            let u = rhs.lsu[i as usize] as i64 * norm_factor + carry;
            carry = u / UNIT as i64;
            buf2[(i + rhs_shift_units) as usize] = (u - carry * UNIT as i64) as i32;
            i += 1;
        }
        debug_assert!(carry == 0, "carry must be 0 in divider normalization");
    }
    let stop = buf2_len - 1; // stop index
    let d0 = buf2[(rhs_non_zero + rhs_shift_units) as usize] as i64; // rhs most significant unit
    let d1 = buf2[(rhs_non_zero + rhs_shift_units - 1) as usize] as i64; // rhs second significant unit
    let mut i = buf1_len + dividend_shift - 1;
    while i >= stop {
        // D3. make the guess on u1
        let u0 = buf1[(i + 1) as usize] as i64;
        let mut u1 = 0;
        if i >= 0 {
            u1 = buf1[i as usize] as i64;
        }
        let u = u0 * UNIT as i64 + u1;
        let mut qhat = u / d0;
        let mut rhat = u - qhat * d0 as i64;
        debug_assert!(qhat < UNIT as i64, "qhat must be less than Unit");
        let mut u2 = 0;
        if i > 0 {
            u2 = buf1[(i - 1) as usize] as i64;
        }
        while qhat * d1 > rhat * UNIT as i64 + u2 {
            // check if qhat can satisfy next unit
            qhat -= 1;
            rhat += d0;
        }
        // D4. multiply and subtract
        let mut carry: i64 = 0;
        let mut borrow: i32 = 0;
        for (b1, b2) in buf1[(i - buf2_len + 1) as usize..(i + 1) as usize]
            .iter_mut()
            .zip(buf2[..buf2_len as usize].iter())
        {
            let mul_v = qhat * *b2 as i64 + carry; // mul
            carry = mul_v / UNIT as i64; // update carry
            let mul_v0 = mul_v - carry * UNIT as i64; // in current unit
            let (sub_v, b0) = sub_with_borrow(*b1, mul_v0 as i32, borrow);
            borrow = b0; // update borrow
            *b1 = sub_v; // update buf1
        }

        borrow += buf1[(i + 1) as usize] - carry as i32;
        if borrow == -1 {
            // qhat is larger, cannot satisfy the whole decimal
            // D6. add back (reverse subtract)
            // should decrease qhat, but not used in remainder
            borrow = 0; // reset borrow to zero
            for b in &mut buf1[(i - buf2_len + 1) as usize..(i + 1) as usize] {
                let (r0, b0) = sub_with_borrow(0, *b, borrow);
                borrow = b0;
                *b = r0;
            }
        }
        buf1[(i + 1) as usize] = 0;
        i -= 1;
    }
    // now we have remainder in buf1
    debug_assert_eq!(buf1[buf1_len as usize], 0, "value must be zero");

    // divide by norm factor, put quotient in result
    let mut rem: i64 = 0;
    let mut res_non_zero: isize = -1;
    i = buf1_len;
    while i > 0 {
        i -= 1;
        let u = rem * UNIT as i64 + buf1[i as usize] as i64;
        if u == 0 {
            continue;
        }
        let q = u / norm_factor;
        rem = u - q * norm_factor; // update remainder
        if q > 0 {
            res.lsu[i as usize] = q as i32; // update result
            if res_non_zero < 0 {
                res_non_zero = i;
            }
        }
    }
    // because we multiply lhs and rhs with identical normFactor, the remainder must be zero.
    debug_assert!(rem == 0, "remainder must be zero");

    if res_non_zero >= rem_frac_units as isize {
        let res_intg_units = (res_non_zero + 1 - rem_frac_units as isize) as i8;
        res.intg = (res_intg_units * DIGITS_PER_UNIT as i8)
            - unit_leading_zeroes(res.lsu[(res_intg_units + rem_frac_units as i8 - 1) as usize])
                as i8;
    } else {
        res.intg = 0;
    }
    res.frac = rem_frac; // keep fractional precision like subtraction
    Ok(())
}

// NOTE: carry can only be 0 or 1
#[inline]
pub(crate) fn add_with_carry(a: i32, b: i32, carry: i32) -> (i32, i32) {
    let r = a + b + carry;
    if r >= UNIT {
        return (r - UNIT, 1);
    }
    (r, 0)
}

// NOTE: borrow can only be 0 or -1
#[inline]
pub(crate) fn sub_with_borrow(a: i32, b: i32, borrow: i32) -> (i32, i32) {
    let r = a - b + borrow;
    if r < 0 {
        return (UNIT + r, -1);
    }
    (r, 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let mut fd1 = FixedDecimal::zero();
        let mut fd2 = FixedDecimal::zero();
        let mut fd3 = FixedDecimal::zero();
        for (input1, input2, expected) in vec![
            ("0", "0", "0"),
            ("0", "1", "1"),
            ("1", "1", "2"),
            ("1", "-1", "0"),
            ("-1", "1", "0"),
            ("-1", "-100", "-101"),
            ("5", "5", "10"),
            ("1.0", "0", "1.0"),
            ("1.0", "0.0", "1.0"),
            ("-1.0", "0.01", "-0.99"),
            ("-0.3", "1.27", "0.97"),
            ("-0.3", "0.2", "-0.1"),
            ("-0.01", "0.001", "-0.009"),
            ("-123", "0.1", "-122.9"),
            ("1", "-12.5", "-11.5"),
            ("-5.0", "5.0", "0.0"),
            ("1.0", "0.1", "1.1"),
            ("1.01", "0.1", "1.11"),
            ("1.00000000001", "1000.01", "1001.01000000001"),
            ("1.234567890", "0.0000000001", "1.2345678901"),
            ("-1.234567890", "0.0000000001", "-1.2345678899"),
            ("10000000000", "1", "10000000001"),
            ("1", "10000000000", "10000000001"),
            ("999999999", "1", "1000000000"),
            ("1", "999999999", "1000000000"),
            ("999999999999999999", "1", "1000000000000000000"),
            ("1", "999999999999999999", "1000000000000000000"),
        ] {
            assert!(fd1.from_ascii_str(input1, true).is_ok());
            assert!(fd2.from_ascii_str(input2, true).is_ok());
            assert!(FixedDecimal::add_to(&fd1, &fd2, &mut fd3).is_ok());
            let actual = fd3.to_string(-1);
            println!("fd3={:?}", fd3);
            println!(
                "{} + {} = {} , intg={}, frac={}",
                fd1.to_string(-1),
                fd2.to_string(-1),
                actual,
                fd3.intg(),
                fd3.frac()
            );
            assert_eq!(&actual[..], expected);
            let actual = (fd1.clone() + fd2.clone()).to_string(-1);
            assert_eq!(&actual[..], expected);
        }
    }

    #[test]
    fn test_sub() {
        let mut fd1 = FixedDecimal::zero();
        let mut fd2 = FixedDecimal::zero();
        let mut fd3 = FixedDecimal::zero();
        for (input1, input2, expected) in vec![
            ("0", "0", "0"),
            ("0", "1", "-1"),
            ("0", "-1", "1"),
            ("-1", "0", "-1"),
            ("1", "0", "1"),
            ("1", "1", "0"),
            ("1", "2", "-1"),
            ("2", "1", "1"),
            ("1", "-1", "2"),
            ("-1", "1", "-2"),
            ("-1", "-100", "99"),
            ("1.0", "0", "1.0"),
            ("1.0", "0.0", "1.0"),
            ("-1.0", "0.01", "-1.01"),
            ("-0.3", "1.27", "-1.57"),
            ("-0.3", "-0.2", "-0.1"),
            ("-0.3", "0.2", "-0.5"),
            ("-0.01", "0.001", "-0.011"),
            ("-123", "0.1", "-123.1"),
            ("1", "-12.5", "13.5"),
            ("-5.0", "5.0", "-10.0"),
            ("1.0", "0.1", "0.9"),
            ("1.01", "0.1", "0.91"),
            ("1.00000000001", "1000.01", "-999.00999999999"),
            ("1.234567890", "0.0000000001", "1.2345678899"),
            ("-1.234567890", "0.0000000001", "-1.2345678901"),
            ("1000000000", "1", "999999999"),
            ("1", "1000000000", "-999999999"),
        ] {
            assert!(fd1.from_ascii_str(input1, true).is_ok());
            assert!(fd2.from_ascii_str(input2, true).is_ok());
            assert!(FixedDecimal::sub_to(&fd1, &fd2, &mut fd3).is_ok());
            let actual = fd3.to_string(-1);
            println!("fd3={:?}", fd3);
            println!(
                "{} - {} = {} , intg={}, frac={}",
                fd1.to_string(-1),
                fd2.to_string(-1),
                actual,
                fd3.intg(),
                fd3.frac()
            );
            assert_eq!(&actual[..], expected);
            let actual = (fd1.clone() - fd2.clone()).to_string(-1);
            assert_eq!(&actual[..], expected);
        }
    }

    #[test]
    fn test_mul() {
        let mut fd1 = FixedDecimal::zero();
        let mut fd2 = FixedDecimal::zero();
        let mut fd3 = FixedDecimal::zero();
        for (input1, input2, expected) in vec![
            ("0", "0", "0"),
            ("0", "1", "0"),
            ("1", "0", "0"),
            ("1", "1", "1"),
            ("1", "2", "2"),
            ("2", "1", "2"),
            ("1", "-1", "-1"),
            ("-1", "1", "-1"),
            ("-1", "-100", "100"),
            ("1.0", "0", "0"),
            ("1.0", "0.0", "0"),
            ("-1.0", "0.01", "-0.010"),
            ("-0.3", "1.27", "-0.381"),
            ("-0.3", "-0.2", "0.06"),
            ("-0.3", "0.2", "-0.06"),
            ("-0.01", "0.001", "-0.00001"),
            ("-0.10", "0.001", "-0.00010"),
            ("-123", "0.1", "-12.3"),
            ("1", "-12.5", "-12.5"),
            ("-5.0", "5.0", "-25.00"),
            ("1.0", "0.1", "0.10"),
            ("1.01", "0.1", "0.101"),
            ("1.00000000001", "1000.01", "1000.0100000100001"),
            ("1.234567890", "0.0000000001", "0.0000000001234567890"),
            ("-1.234567890", "0.0000000001", "-0.0000000001234567890"),
        ] {
            assert!(fd1.from_ascii_str(input1, true).is_ok());
            assert!(fd2.from_ascii_str(input2, true).is_ok());
            assert!(FixedDecimal::mul_to(&fd1, &fd2, &mut fd3).is_ok());
            let actual = fd3.to_string(-1);
            println!("fd3={:?}", fd3);
            println!(
                "{} * {} = {} , intg={}, frac={}",
                fd1.to_string(-1),
                fd2.to_string(-1),
                actual,
                fd3.intg(),
                fd3.frac()
            );
            assert_eq!(&actual[..], expected);
            let actual = (fd1.clone() * fd2.clone()).to_string(-1);
            assert_eq!(&actual[..], expected);
        }
    }

    #[test]
    fn test_div() {
        let mut fd1 = FixedDecimal::zero();
        let mut fd2 = FixedDecimal::zero();
        let mut fd3 = FixedDecimal::zero();
        for (input1, input2, expected) in vec![
            ("0", "1", "0"),
            ("1", "1", "1.000000000"), // incremental frac is already multiple of 9
            ("1", "1", "1.000000000"),
            ("1", "2", "0.500000000"),
            ("2", "1", "2.000000000"),
            ("1", "-1", "-1.000000000"),
            ("-1", "1", "-1.000000000"),
            ("-1", "-100", "0.010000000"),
            ("100", "1", "100.000000000"),
            ("100", "100", "1.000000000"),
            ("0.000000002", "1", "0.000000002000000000"),
            ("1.0", "2", "0.500000000"),
            ("1.0", "2.0", "0.500000000000000000"), // two fractional units
            ("-1", "0.01", "-100.000000000"),
            ("0.27", "0.3", "0.900000000000000000"),
            ("-0.3", "-0.2", "1.500000000000000000"),
            ("0.3", "0.7", "0.428571428571428571"),
            ("0.6", "0.9", "0.666666666666666666"),
            ("-0.3", "0.2", "-1.500000000000000000"),
            ("1000000000.1", "7", "142857142.871428571"),
            ("1000000000.1", "9", "111111111.122222222"),
            ("101000000000.1", "7", "14428571428.585714285"),
            ("101000000000.1", "7.1", "14225352112.690140845070422535"),
            ("101000000000.1", "5", "20200000000.020000000"),
            ("101000000000.1", "5.0", "20200000000.020000000000000000"),
            ("100.10000000001", "7", "14.300000000001428571"),
            ("100.10000000001", "7.0", "14.300000000001428571428571428"),
            ("100.1", "7.0000000001", "14.299999999795714285717204081"),
            ("205.6", "9.5000000001", "21.642105262930083102495472809"),
            (
                "2000000005.1",
                "7.5000000001",
                "266666667.343111111102091851851972108",
            ),
            ("1.2", "0.7", "1.714285714285714285"),
            ("1.22", "0.77", "1.584415584415584415"),
            ("1.222", "0.777", "1.572715572715572715"),
            ("1.2222", "0.7777", "1.571557155715571557"),
            ("1.22222", "0.77777", "1.571441428700001285"),
            ("1.222222", "0.777777", "1.571429857144142858"),
            ("1.2222222", "0.7777777", "1.571428700000012857"),
            ("1.22222222", "0.77777777", "1.571428584285714414285715571"),
            (
                "1.222222222",
                "0.777777777",
                "1.571428572714285715571428572",
            ),
            ("9.8", "1", "9.800000000"),
            ("98.7", "1.2", "82.250000000000000000"),
            ("987.6", "12.3", "80.292682926829268292"),
            ("9876.5", "123.4", "80.036466774716369529"),
            ("98765.4", "1234.5", "80.004374240583232077"),
            ("987654.3", "12345.6", "80.000510303265940902"),
            ("9876543.2", "123456.7", "80.000058320042573631"),
            ("98765432.1", "1234567.8", "80.000006561000538002"),
            ("987654321.1", "12345678.9", "80.000000737100006707"),
            ("987654321.12", "12345678.99", "80.000000155520000281"),
            ("987654321.123", "12345678.998", "80.000000103923000120"),
            ("987654321.1234", "12345678.9987", "80.000000099419400109"),
            ("987654321.12345", "12345678.99876", "80.000000099034650108"),
            (
                "987654321.123456",
                "12345678.998765",
                "80.000000099002736108",
            ),
            (
                "987654321.1234567",
                "12345678.9987654",
                "80.000000099000200808",
            ),
            (
                "987654321.12345678",
                "12345678.99876543",
                "80.000000099000012888900031007",
            ),
            (
                "987654321.123456789",
                "12345678.998765432",
                "80.000000099000000657900001515",
            ),
            ("-9.8", "1", "-9.800000000"),
            ("-98.7", "1.2", "-82.250000000000000000"),
            ("-987.6", "12.3", "-80.292682926829268292"),
            ("-9876.5", "123.4", "-80.036466774716369529"),
            ("-98765.4", "1234.5", "-80.004374240583232077"),
            ("-987654.3", "12345.6", "-80.000510303265940902"),
            ("-9876543.2", "123456.7", "-80.000058320042573631"),
            ("-98765432.1", "1234567.8", "-80.000006561000538002"),
            ("-987654321.1", "12345678.9", "-80.000000737100006707"),
            ("-987654321.12", "12345678.99", "-80.000000155520000281"),
            ("-987654321.123", "12345678.998", "-80.000000103923000120"),
            ("-987654321.1234", "12345678.9987", "-80.000000099419400109"),
            (
                "-987654321.12345",
                "12345678.99876",
                "-80.000000099034650108",
            ),
            (
                "-987654321.123456",
                "12345678.998765",
                "-80.000000099002736108",
            ),
            (
                "-987654321.1234567",
                "12345678.9987654",
                "-80.000000099000200808",
            ),
            (
                "-987654321.12345678",
                "12345678.99876543",
                "-80.000000099000012888900031007",
            ),
            (
                "-987654321.123456789",
                "12345678.998765432",
                "-80.000000099000000657900001515",
            ),
            ("0.170511", "-353390023.459963", "-0.000000000482500887"),
            ("0.170511", "-353390023", "-0.000000000482500888"),
            ("0.1", "300000000", "0.000000000"),
            ("0.1", "300000000.0", "0.000000000333333333"),
            ("0.1", "3000000000", "0.000000000"),
            ("0.1", "3000000000.0", "0.000000000033333333"),
            ("0.0000000001", "300000000", "0.000000000000000000"),
            (
                "0.0000000001",
                "300000000.0",
                "0.000000000000000000333333333",
            ),
            ("0.0000000001", "3000000000", "0.000000000000000000"),
            (
                "0.0000000001",
                "3000000000.0",
                "0.000000000000000000033333333",
            ),
            ("1", "300000000", "0.000000003"),
            ("1", "300000000.0", "0.000000003"),
            ("1", "3000000000", "0.000000000"),
            ("1", "3000000000.0", "0.000000000"),
            ("1.0", "300000000", "0.000000003"),
            ("1.0", "300000000.0", "0.000000003333333333"),
            ("1.0", "3000000000", "0.000000000"),
            ("1.0", "3000000000.0", "0.000000000333333333"),
            ("0.4", "0.000000003", "133333333.333333333333333333"),
            (
                "0.4",
                "0.0000000003",
                "1333333333.333333333333333333333333333",
            ),
            ("0.2", "0.000000003", "66666666.666666666666666666"),
            (
                "0.2",
                "0.0000000003",
                "666666666.666666666666666666666666666",
            ),
            ("400000000", "300000000", "1.333333333"),
            ("400000000.0", "300000000.0", "1.333333333333333333"),
            ("4000000000", "3000000000", "1.333333333"),
            ("4000000000.0", "3000000000.0", "1.333333333333333333"),
            ("200000000", "300000000", "0.666666666"),
            ("200000000.0", "300000000.0", "0.666666666666666666"),
            ("2000000000", "3000000000", "0.666666666"),
            ("2000000000.0", "3000000000.0", "0.666666666666666666"),
            (
                "400000000",
                "0.000000003",
                "133333333333333333.333333333333333333",
            ),
            (
                "4000000000",
                "0.000000003",
                "1333333333333333333.333333333333333333",
            ),
            ("1", "500000000.1", "0.000000001"),
        ] {
            assert!(fd1.from_ascii_str(input1, true).is_ok());
            assert!(fd2.from_ascii_str(input2, true).is_ok());
            assert!(FixedDecimal::div_to(&fd1, &fd2, &mut fd3, DIV_INCR_FRAC).is_ok());
            let actual = fd3.to_string(-1);
            println!("fd3={:?}", fd3);
            println!(
                "{} / {} = {} , intg={}, frac={}",
                fd1.to_string(-1),
                fd2.to_string(-1),
                actual,
                fd3.intg(),
                fd3.frac()
            );
            assert_eq!(&actual[..], expected);
            let actual = (fd1.clone() / fd2.clone()).to_string(-1);
            assert_eq!(&actual[..], expected);
        }
    }

    #[test]
    fn test_rem() {
        let mut fd1 = FixedDecimal::zero();
        let mut fd2 = FixedDecimal::zero();
        let mut fd3 = FixedDecimal::zero();
        for (input1, input2, expected) in vec![
            ("0", "1", "0"),
            ("1", "1", "0"),
            ("1", "2", "1"),
            ("2", "1", "0"),
            ("1000000001", "2", "1"),
            ("-1000000001", "2", "-1"),
            ("-1", "2", "-1"),
            ("-1", "-2", "-1"),
            ("1", "-2", "1"),
            ("-1", "-100", "-1"),
            ("100", "3", "1"),
            ("100", "1001", "100"),
            ("0.2", "1", "0.2"),
            ("0.02", "1", "0.02"),
            ("0.002", "1", "0.002"),
            ("0.0002", "1", "0.0002"),
            ("0.00002", "1", "0.00002"),
            ("0.000002", "1", "0.000002"),
            ("0.0000002", "1", "0.0000002"),
            ("0.00000002", "1", "0.00000002"),
            ("0.000000002", "1", "0.000000002"),
            ("0.2", "1.0", "0.2"),
            ("0.2", "1.00", "0.20"),
            ("0.2", "1.000", "0.200"),
            ("0.2", "1.0000", "0.2000"),
            ("0.2", "1.00000", "0.20000"),
            ("0.2", "1.000000", "0.200000"),
            ("0.2", "1.0000000", "0.2000000"),
            ("0.2", "1.00000000", "0.20000000"),
            ("0.2", "1.000000000", "0.200000000"),
            ("-0.2", "1.0", "-0.2"),
            ("-0.2", "1.00", "-0.20"),
            ("-0.2", "1.000", "-0.200"),
            ("-0.2", "1.0000", "-0.2000"),
            ("-0.2", "1.00000", "-0.20000"),
            ("-0.2", "1.000000", "-0.200000"),
            ("-0.2", "1.0000000", "-0.2000000"),
            ("-0.2", "1.00000000", "-0.20000000"),
            ("-0.2", "1.000000000", "-0.200000000"),
            ("-0.3", "-0.2", "-0.1"),
            ("0.3", "0.2", "0.1"),
            ("0.3", "-0.2", "0.1"),
            ("-0.3", "0.2", "-0.1"),
            ("-0.3", "-0.7", "-0.3"),
            ("0.3", "-0.7", "0.3"),
            ("-0.3", "0.7", "-0.3"),
            ("0.3", "0.7", "0.3"),
            ("1000000000.1", "7", "6.1"),
            ("1000000000.1", "70", "20.1"),
            ("101000000000.1", "2000000000", "1000000000.1"),
            ("1000000000.1", "9", "1.1"),
            ("1000000000.1", "9.00", "1.10"),
            ("100.10000000001", "7", "2.10000000001"),
            ("101000000000.1", "7.1", "4.9"),
            ("101000000000.1", "5", "0.1"),
            ("101000000000.1", "5.291", "0.201"),
            ("100.1", "7.0000000001", "2.0999999986"),
            ("205.6", "9.5000000001", "6.0999999979"),
            ("2000000005.1", "7.5000000001", "2.5733333333"),
            ("1.2", "0.7", "0.5"),
            ("1.22", "0.77", "0.45"),
            ("1.222", "0.777", "0.445"),
            ("1.2222", "0.7777", "0.4445"),
            ("1.22222", "0.77777", "0.44445"),
            ("1.222222", "0.777777", "0.444445"),
            ("1.2222222", "0.7777777", "0.4444445"),
            ("1.22222222", "0.77777777", "0.44444445"),
            ("1.222222222", "0.777777777", "0.444444445"),
            ("9.8", "1", "0.8"),
            ("98.7", "1.2", "0.3"),
            ("987.6", "1.23", "1.14"),
            ("9876.5", "1.234", "0.798"),
            ("98765.4", "1.2345", "0.4620"),
            ("987654.3", "1.23456", "0.12720"),
            ("9876543.2", "1.234567", "1.027165"),
            ("98765432.1", "1.2345678", "0.6925932"),
            ("987654321.1", "1.23456789", "0.45802477"),
            ("987654321.12", "1.234567899", "0.685432101"),
            ("987654321.123", "1.2345678998", "0.0484321002"),
            ("987654321.1234", "1.23456789987", "1.22740000000"),
            ("987654321.12345", "1.234567899876", "1.222650000000"),
            ("987654321.123456", "1.2345678998765", "1.2222560000000"),
            ("987654321.1234567", "1.23456789987654", "1.22222470000000"),
            (
                "987654321.12345678",
                "1.234567899876543",
                "1.222222380000000",
            ),
            (
                "987654321.123456789",
                "1.2345678998765432",
                "1.2222222290000000",
            ),
            ("-9.8", "1", "-0.8"),
            ("-98.7", "1.2", "-0.3"),
            ("-987.6", "1.23", "-1.14"),
            ("-9876.5", "1.234", "-0.798"),
            ("-98765.4", "1.2345", "-0.4620"),
            ("-987654.3", "1.23456", "-0.12720"),
            ("-9876543.2", "1.234567", "-1.027165"),
            ("-98765432.1", "1.2345678", "-0.6925932"),
            ("-987654321.1", "1.23456789", "-0.45802477"),
            ("-987654321.12", "1.234567899", "-0.685432101"),
            ("-987654321.123", "1.2345678998", "-0.0484321002"),
            ("-987654321.1234", "1.23456789987", "-1.22740000000"),
            ("-987654321.12345", "1.234567899876", "-1.222650000000"),
            ("-987654321.123456", "1.2345678998765", "-1.2222560000000"),
            (
                "-987654321.1234567",
                "1.23456789987654",
                "-1.22222470000000",
            ),
            (
                "-987654321.12345678",
                "1.234567899876543",
                "-1.222222380000000",
            ),
            (
                "-987654321.123456789",
                "1.2345678998765432",
                "-1.2222222290000000",
            ),
            ("0.170511", "-353390023.459963", "0.170511"),
            ("-353390023.459963", "0.170511", "-0.060946"),
            ("0.170511", "-353390023", "0.170511"),
            ("-353390023", "0.170511", "-0.112516"),
            ("0.4", "0.000000003", "0.000000001"),
            ("0.4", "0.0000000003", "0.0000000001"),
            ("0.2", "0.000000003", "0.000000002"),
            ("0.2", "0.0000000003", "0.0000000002"),
            ("1000000000000000001", "70298007", "68215565"),
            ("1000000000000000001", "0.70298007", "0.07924142"),
            ("1000000000000000001", "500000000.1", "300000001.1"),
            ("0.1", "0.20000000001", "0.10000000000"),
        ] {
            assert!(fd1.from_ascii_str(input1, true).is_ok());
            assert!(fd2.from_ascii_str(input2, true).is_ok());
            assert!(FixedDecimal::rem_to(&fd1, &fd2, &mut fd3).is_ok());
            let actual = fd3.to_string(-1);
            println!("fd3={:?}", fd3);
            println!(
                "{} % {} = {} , intg={}, frac={}",
                fd1.to_string(-1),
                fd2.to_string(-1),
                actual,
                fd3.intg(),
                fd3.frac()
            );
            assert_eq!(&actual[..], expected);
            let actual = (fd1.clone() % fd2.clone()).to_string(-1);
            assert_eq!(&actual[..], expected);
        }
    }
}
