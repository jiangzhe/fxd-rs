use crate::decimal::{FixedDecimal, MAX_DIGITS, MAX_FRAC, DIGITS_PER_UNIT};
use crate::error::{Result, Error};
use crate::utils::{BIN_HIGH_UNIT, BIN_MID_UNIT, POW10, BIN2CHAR, d3_to_str, get_units, mod9};

use std::cmp::Ordering;

impl FixedDecimal {
    
    /// parses numeric string and updates decimal value
    #[inline]
    pub fn from_ascii_str(&mut self, s: &str, reset: bool) -> Result<()> {
        self.from_bytes_str(s.as_bytes(), reset)
    }

    /// parses numeric string and updates decimal value
    pub fn from_bytes_str(&mut self, bs: &[u8], reset: bool) -> Result<()> {
        if reset {
            self.set_zero();
        }
        if bs.is_empty() {
            return Err(Error::ConversionSyntax)
        }
        let mut exp = 0; // exponent
        let mut d = 0; // total digits of input
        let mut dotchar = -1; // position of dot character
        let mut last = -1; // position of last digit
        let mut cfirst = 0; // position of first digit
        let mut neg = false; // if input is negative
        let mut more_to_proc = false;
        let mut i: isize = -1;
        for c in bs.iter() {
            i += 1;
            if *c >= b'0' && *c <= b'9' { // test for digit char
                last = i as isize;
                d += 1;
                continue;
            }
            if *c == b'.' && dotchar == -1 {
                dotchar = i as isize;
                if i == cfirst {
                    cfirst += 1; // first digit must follow
                }
                continue;
            }
            if i == 0 { // first in string
                if *c == b'-' { // valid - sign
                    cfirst += 1;
                    neg = true;
                    continue
                }
                if *c == b'+' { // valid + sign
                    cfirst += 1;
                    continue
                }
            }
            // c is not a digit, or a valid '+', '-' or '.'
            more_to_proc = true; // indicate more data to process
            break
        }
        if last == -1 { // no digits yet
            // decided not to support Infinities and NaNs now
            return Err(Error::ConversionSyntax);
        }
        if more_to_proc {
            // has some digits; exponent is only valid sequence now
            let mut nege = false; // negative exponent
            let mut c = bs[i as usize];
            if c != b'e' && c != b'E' {
                return Err(Error::ConversionSyntax);
            }
            // found 'e' or 'E'
            // sign no longer required
            i += 1;
            if i as usize == bs.len() { // to (possible) sign
                return Err(Error::ConversionSyntax);
            }
            c = bs[i as usize];
            if c == b'-' {
                nege = true;
                i += 1;
                if i as usize == bs.len() {
                    return Err(Error::ConversionSyntax)
                }
            } else if c == b'+' {
                i += 1;
                if i as usize == bs.len() {
                    return Err(Error::ConversionSyntax)
                }
            }
            loop {
                c = bs[i as usize];
                if c == b'0' && i as usize != bs.len()-1 { // strip insignificant zeroes
                    i += 1;
                } else {
                    break;
                }
            }
            let firstexp = i; // save exponent digit place
            for c in bs[i as usize..].iter() {
                if *c < b'0' || *c > b'9' { // not a digit
                    return Err(Error::ConversionSyntax);
                }
                exp = exp*10 + *c as i32 - '0' as i32;
                i += 1
            }
            // maximum exponent is 65, with sign at most 4 chars
            if i >= firstexp + 4 {
                return Err(Error::ConversionSyntax);
            }
            if (!nege && exp as usize > MAX_DIGITS) || (nege && exp as usize > MAX_FRAC) {
                return Err(Error::ConversionSyntax);
            }
            if nege {
                exp = -exp;
            }
        }
        // Here when whole string has been inspected; syntax is good
        // cfirst->first digit(never dot), last->last digit(ditto)
        let frac_digits = if dotchar == -1 || last < dotchar { // no dot found or dot is last char
            0
        } else {
            last - dotchar
        };
        // apply exponent to frac
        let (frac, digits, mut heading_zeroes) = {
            let mut frac = frac_digits as i8 - exp as i8;
            let digits;
            let mut heading_zeroes = 0i8;
            match frac.cmp(&0) {
                Ordering::Equal => digits = d,
                Ordering::Greater => {
                    if d > frac { // have both integral and fractional part
                        digits = d;
                    } else { // only fractional part
                        digits = frac;
                        heading_zeroes = frac - d;
                    }
                }
                Ordering::Less => {
                    digits = d - frac;
                    frac = 0;
                }
            }
            (frac, digits, heading_zeroes)
        };

        if digits as usize > MAX_DIGITS {
            return Err(Error::ConversionSyntax);
        }

        // units of integral part and fractional part
        let intg_units = get_units(digits - frac);
        let frac_units = get_units(frac);
        let mut up = (intg_units + frac_units -1) as isize; // lsu unit index, from highest to lowest
        // reset i as parse index
        i = cfirst;
        if intg_units > 0 {
            let mut out = 0;
            let mut cut = digits -frac;
            let mut c;
            loop {
                c = bs[i as usize];
                if c == b'.' { // ignore '.', this may be caused by exponent normalization
                    i += 1;
                    continue;
                }
                out = out*10 + c as i32 - b'0' as i32;
                cut -= 1;
                if cut == 0 {
                    break // nothing for this unit
                }
                if i == last { // no more digits to read, adjust out if cut > 0
                    break
                }
                if mod9(cut) > 0 {
                    i += 1;
                    continue;
                }
                self.lsu[up as usize] = out; // write out
                up -= 1; // prepare for unit below
                out = 0;
                i += 1;
            }
            // input integral digits processed.
            // increment i for frac processing
            i += 1;
            // handle cut > 0, e.g. "1E2", "1E20"
            let re = mod9(cut);
            self.lsu[up as usize] = out * POW10[re];
            up -= 1;
            cut -= re as i8;
            while cut > 0 {
                self.lsu[up as usize] = 0;
                up -= 1;
                cut -= DIGITS_PER_UNIT as i8;
            }
        }
        if frac_units > 0 {
            let mut out = 0; // accumulator in unit
            let mut cut = DIGITS_PER_UNIT as i8;
            while heading_zeroes >= DIGITS_PER_UNIT as i8 { // current unit filled all zeroes
                self.lsu[up as usize] = 0;
                up -= 1;
                heading_zeroes -= DIGITS_PER_UNIT as i8;
            }
            cut -= heading_zeroes as i8;
            // parse fractional part
            let mut c;
            loop {
                c = bs[i as usize];
                if c == b'.' { // ignore '.', this may be caused by exponent normalization
                    i += 1;    
                    continue;
                }
                cut -= 1;
                out += (c as i32 - '0' as i32) * POW10[cut as usize];
                if i == last { // done
                    break;
                }
                if cut > 0 {
                    i += 1;
                    continue; // more for this unit
                }
                self.lsu[up as usize] = out; // write unit
                up -= 1; // prepare for unit below
                cut = DIGITS_PER_UNIT as i8;
                out = 0;
                i += 1;
            }
            self.lsu[up as usize] = out; // write lsu
        }
        self.intg = digits - frac;
        self.frac = frac;
        if neg {
            self.set_neg();
        }
        Ok(())
    }

    /// converts this decimal to string format.
    /// frac specifies the fractional precision of the output string.
    /// if frac < 0, will output all fractional digits.
    /// if frac >= 0, will truncate to given digits.
    pub fn to_string(&self, frac: isize) -> String {
        let mut bs = Vec::new();
        self.append_str_buf(&mut bs, frac);
        unsafe { String::from_utf8_unchecked(bs) }
    }

    /// appends formatted string to given buffer.
    /// if frac < 0, will output all fractional digits.
    /// if frac >= 0, will truncate to given digits.
    pub fn append_str_buf(&self, buf: &mut Vec<u8>, frac: isize) {
        if self.is_neg() {
            buf.push(b'-');
        }
        let intg_units = self.intg_units();
        let frac_units = self.frac_units();
        let mut up = (intg_units + frac_units) as isize - 1;
        if intg_units > 0 { // intg part
            let mut elim_heading_zeroes = true;
            while up >= frac_units as isize {
                // for each unit, we divide into 3 3-digit parts and print them each time
                let mut u = self.lsu[up as usize];
                up -= 1;
                if elim_heading_zeroes && u == 0 { // rare case, integral part exists but zero
                    continue;
                }
                // XXX,xxx,xxx
                if u >= BIN_HIGH_UNIT {
                    let high = u / BIN_HIGH_UNIT;
                    elim_heading_zeroes = d3_to_str(high, buf, elim_heading_zeroes);
                    u -= high * BIN_HIGH_UNIT;
                } else {
                    elim_heading_zeroes = d3_to_str(0, buf, elim_heading_zeroes);
                }
                // xxx,XXX,xxx
                if u >= BIN_MID_UNIT {
                    let mid = u / BIN_MID_UNIT;
                    elim_heading_zeroes = d3_to_str(mid, buf, elim_heading_zeroes);
                    u -= mid * BIN_MID_UNIT;
                } else {
                    elim_heading_zeroes = d3_to_str(0, buf, elim_heading_zeroes);
                }
                // xxx,xxx,XXX
                elim_heading_zeroes = d3_to_str(u, buf, elim_heading_zeroes);
            }
            if elim_heading_zeroes { // integral units are all zeroes
                buf.push(b'0');
            }
        } else {
            buf.push(b'0'); // append 0 if no integral part
        }
        if frac == 0 { // no fractional digits
            return;
        }
        if frac < 0 { // frac is not specified
            if frac_units > 0 { // fractional part exists
                buf.push(b'.');
                let frac_digits = self.frac();
                self.append_frac(buf, up as isize, frac_digits);
            }
            return;
        }
        // frac is specified
        if frac_units == 0 { // no fractional part
            buf.push(b'.');
            for _ in 0..frac {
                buf.push(b'0');
            }
            return;
        }
        let full_frac_digits = (frac_units * DIGITS_PER_UNIT) as isize;
        if frac <= full_frac_digits { // this decimal supports specified frac
            buf.push(b'.');
            self.append_frac(buf, up, frac as i8);
            return;
        }
        // this decimal does not have sufficient precision
        buf.push(b'.');
        self.append_frac(buf, up, full_frac_digits as i8);
        // add extra zeroes
        for _ in 0..frac-full_frac_digits {
            buf.push(b'0');
        }
    }

    fn append_frac(&self, buf: &mut Vec<u8>, mut up: isize, mut frac_digits: i8) {
        while up >= 0 {
            // for each unit, divide into 3 3-digit parts and print.
            // need to check whether the digits exceeds fractional length.
            if frac_digits == 0 { // already reached the fractional precision limit
                return
            }
            let mut u = self.lsu[up as usize];
            up -= 1;
            if u == 0 { // all zeros
                if frac_digits > DIGITS_PER_UNIT as i8 { // digit number larger than 
                    buf.extend_from_slice(b"000000000");
                    frac_digits -= DIGITS_PER_UNIT as i8;
                    continue
                }
                // append n zeroes and return
                for _ in 0..frac_digits { 
                    buf.push(b'0');
                }
                return
            }
            // non-zero
            // XXX,xxx,xxx
            if u >= BIN_HIGH_UNIT {
                let high = u / BIN_HIGH_UNIT;
                let start_idx = high as usize * 4 + 1;
                if frac_digits < 3 {
                    buf.extend_from_slice(&BIN2CHAR[start_idx..start_idx+frac_digits as usize]);
                    return;
                }
                buf.extend_from_slice(&BIN2CHAR[start_idx..start_idx+3]);
                frac_digits -= 3;
                if frac_digits == 0 {
                    return;
                }
                u -= high * BIN_HIGH_UNIT;
            } else {
                if frac_digits < 3 {
                    for _ in 0..frac_digits {
                        buf.push(b'0');
                    }
                    return;
                }
                buf.extend_from_slice(b"000");
                frac_digits -= 3;
                if frac_digits == 0 {
                    return;
                }
            }
            // xxx,XXX,xxx
            if u >= BIN_MID_UNIT {
                let mid = u / BIN_MID_UNIT;
                let start_idx = mid as usize * 4 + 1;
                if frac_digits < 3 {
                    buf.extend_from_slice(&BIN2CHAR[start_idx..start_idx+frac_digits as usize]);
                    return;
                }
                buf.extend_from_slice(&BIN2CHAR[start_idx..start_idx+3]);
                frac_digits -= 3;
                if frac_digits == 0 {
                    return;
                }
                u -= mid * BIN_MID_UNIT;
            } else {
                if frac_digits < 3 {
                    for _ in 0..frac_digits {
                        buf.push(b'0');
                    }
                    return;
                }
                buf.extend_from_slice(b"000");
                frac_digits -= 3;
                if frac_digits == 0 {
                    return;
                }
            }
            // xxx,xxx,XXX
            let start_idx = u as usize * 4 + 1;
            if frac_digits < 3 {
                buf.extend_from_slice(&BIN2CHAR[start_idx..start_idx+frac_digits as usize]);
                return;
            }
            buf.extend_from_slice(&BIN2CHAR[start_idx..start_idx+3]);
            frac_digits -= 3;
            if frac_digits == 0 {
                return;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    
    use super::*;

    #[test]
    fn test_from_str() {
        let mut fd = FixedDecimal::zero();
        for s in vec![
            "0",
            "1",
            "-1",
            "123",
            "123456789012345",
            "0.1",
            "0.123",
            "1.0",
            "-1.0",
            "1E1",
            "1E+2",
            "1.0E-2",
            "1.0e10",
            "1.0e02",
        ] {
            assert!(fd.from_ascii_str(s, true).is_ok());
            println!("result={:?}", fd);
        }
    }

    #[test]
    fn test_from_str_err() {
        let mut fd = FixedDecimal::zero();
        for s in vec![
            "",
            "abc",
            ".",
            ".a",
            ".NaN",
            "N",
            "Nb",
            "Na",
            "Nab",
            "NaNx",
            "0x",
            "0E",
            ".1e+",
            ".1e-",
            ".1e+f",
            ".1e12345",
            ".1e-200",
            "1234567890123456789012345678901234567890123456789012345678901234567890",
        ] {
            assert!(fd.from_ascii_str(s, true).is_err());
        }
    }

    #[test]
    fn test_to_string() {
        let mut fd = FixedDecimal::zero();
        for (input, frac, expected) in vec![
                ("0", -1, "0"),
                ("1", -1, "1"),
                ("-1", -1, "-1"),
                ("+1", -1, "1"),
                ("123", -1, "123"),
                ("123456789012345", -1, "123456789012345"),
                ("0.0", -1, "0.0"),
                ("0.100", -1, "0.100"),
                ("0.1", -1, "0.1"),
                ("0.12345678901234567890", -1, "0.12345678901234567890"),
                ("0.123", -1, "0.123"),
                ("1.0", -1, "1.0"),
                ("-1.0", -1, "-1.0"),
                ("1e0", -1, "1"),
                ("1E1", -1, "10"),
                ("1E+2", -1, "100"),
                ("1.0E-2", -1, "0.010"),
                ("1.0e10", -1, "10000000000"),
                ("1.2345e20", -1, "123450000000000000000"),
                ("5.4433e4", -1, "54433"),
                ("5.4433e3", -1, "5443.3"),
                ("5.4433e2", -1, "544.33"),
                ("5.4433e1", -1, "54.433"),
                ("5.4433e0", -1, "5.4433"),
                ("5.4433e-1", -1, "0.54433"),
                ("5.4433e-2", -1, "0.054433"),
                ("5.4433e-3", -1, "0.0054433"),
                ("5.4433e-5", -1, "0.000054433"),
                ("5.4433e-6", -1, "0.0000054433"),
                ("5.4433e-7", -1, "0.00000054433"),
                ("5.4433e-8", -1, "0.000000054433"),
                ("5.4433e-9", -1, "0.0000000054433"),
                ("5.4433e-10", -1, "0.00000000054433"),
                ("5.4433e-11", -1, "0.000000000054433"),
                ("5.4433e-20", -1, "0.000000000000000000054433"),
                ("123456789.123456789", 0, "123456789"),
                ("123456789.123456789", 1, "123456789.1"),
                ("123456789.123456789", 2, "123456789.12"),
                ("123456789.123456789", 3, "123456789.123"),
                ("123456789.123456789", 4, "123456789.1234"),
                ("123456789.123456789", 5, "123456789.12345"),
                ("123456789.123456789", 6, "123456789.123456"),
                ("123456789.123456789", 7, "123456789.1234567"),
                ("123456789.123456789", 8, "123456789.12345678"),
                ("123456789.123456789", 9, "123456789.123456789"),
                ("123456789.123456789", 10, "123456789.1234567890"),
                ("1.234000567", 4, "1.2340"),
                ("1.234000567", 5, "1.23400"),
                ("1.234000567", 6, "1.234000"),
                ("1.021", 2, "1.02"),
                ("1.00021", 2, "1.00"),
                ("123", 2, "123.00"),
        ] {
            assert!(fd.from_ascii_str(&input, true).is_ok());
            println!("fd={:?}", fd);
            let actual = fd.to_string(frac);
            assert_eq!(expected, &actual[..]);
        }
    }

}