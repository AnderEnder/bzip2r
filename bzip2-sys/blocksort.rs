use crate::private_ffi::{fallbackSort, mainSort, EState, BZ_N_OVERSHOOT};
use compress::{assertd, asserth};
use std::slice::from_raw_parts_mut;

#[no_mangle]
pub extern "C" fn fallbackSimpleSort(fmap: *mut u32, eclass: *mut u32, lo: i32, hi: i32) {
    let fmap = unsafe { from_raw_parts_mut(fmap, (hi + 1) as usize) };
    let max_index = fmap
        .iter()
        .skip(lo as usize)
        .max()
        .cloned()
        .unwrap_or_default() as usize;
    let eclass = unsafe { from_raw_parts_mut(eclass, max_index + 1) };
    fallback_simple_sort(fmap, eclass, lo, hi);
}

// some overflow is present
fn fallback_simple_sort(fmap: &mut [u32], eclass: &mut [u32], lo: i32, hi: i32) {
    if lo == hi {
        return;
    }

    if hi - lo > 3 {
        let mut i = hi - 4;
        while i >= lo {
            let tmp = fmap[i as usize];
            let ec_tmp = eclass[tmp as usize];

            let mut j = i + 4;
            while j <= hi && ec_tmp > eclass[fmap[j as usize] as usize] {
                fmap[(j - 4) as usize] = fmap[j as usize];
                j += 4
            }
            fmap[(j - 4) as usize] = tmp;
            i -= 1;
        }
    }

    let mut i = hi - 1;
    while i >= lo {
        let tmp = fmap[i as usize];
        let ec_tmp = eclass[tmp as usize];
        let mut j = i + 1;
        while j <= hi && ec_tmp > eclass[fmap[j as usize] as usize] {
            fmap[(j - 1) as usize] = fmap[j as usize];
            j += 1;
        }
        fmap[(j - 1) as usize] = tmp;
        i -= 1;
    }
}

#[no_mangle]
pub extern "C" fn mmed3(mut a: u8, mut b: u8, c: u8) -> u8 {
    if a > b {
        let t = a;
        a = b;
        b = t;
    };
    if b > c {
        b = c;
        if a > b {
            b = a;
        }
    }
    return b;
}

#[no_mangle]
pub extern "C" fn mainGtU(
    mut i1: u32,
    mut i2: u32,
    block: *mut u8,
    quadrant: *mut u16,
    nblock: u32,
    budget: *mut i32,
) -> u8 {
    let block = unsafe { from_raw_parts_mut(block, (nblock + 8) as usize) };
    let quadrant = unsafe { from_raw_parts_mut(quadrant, (nblock + 8) as usize) };
    main_gtu(i1, i2, block, quadrant, nblock, budget) as u8
}

fn main_gtu(
    mut i1: u32,
    mut i2: u32,
    block: &mut [u8],
    quadrant: &mut [u16],
    nblock: u32,
    budget: *mut i32,
) -> bool {
    assertd(i1 != i2, "mainGtU");
    /* 1 */
    let mut c1 = block[i1 as usize];
    let mut c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 2 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 3 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 4 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 5 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 6 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 7 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 8 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 9 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 10 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 11 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;
    /* 12 */
    c1 = block[i1 as usize];
    c2 = block[i2 as usize];
    if c1 != c2 {
        return c1 > c2;
    }
    i1 += 1;
    i2 += 1;

    let mut k = (nblock as i32) - 8;

    loop {
        /* 1 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }

        let mut s1 = quadrant[i1 as usize];
        let mut s2 = quadrant[i2 as usize];
        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;
        /* 2 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }
        s1 = quadrant[i1 as usize];
        s2 = quadrant[i2 as usize];
        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;
        /* 3 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }
        s1 = quadrant[i1 as usize];
        s2 = quadrant[i2 as usize];
        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;
        /* 4 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }
        s1 = quadrant[i1 as usize];
        s2 = quadrant[i2 as usize];
        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;
        /* 5 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }
        s1 = quadrant[i1 as usize];
        s2 = quadrant[i2 as usize];
        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;
        /* 6 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }
        s1 = quadrant[i1 as usize];
        s2 = quadrant[i2 as usize];
        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;
        /* 7 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }
        s1 = quadrant[i1 as usize];
        s2 = quadrant[i2 as usize];
        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;
        /* 8 */
        c1 = block[i1 as usize];
        c2 = block[i2 as usize];
        if c1 != c2 {
            return c1 > c2;
        }
        s1 = quadrant[i1 as usize];
        s2 = quadrant[i2 as usize];

        if s1 != s2 {
            return s1 > s2;
        }
        i1 += 1;
        i2 += 1;

        if i1 >= nblock {
            i1 -= nblock
        }

        if i2 >= nblock {
            i2 -= nblock
        }

        k -= 8;
        unsafe {
            *budget -= 1;
        }

        if k < 0 {
            break;
        }
    }

    return false;
}

const incs: [i32; 14] = [
    1, 4, 13, 40, 121, 364, 1093, 3280, 9841, 29524, 88573, 265720, 797161, 2391484,
];

#[no_mangle]
pub extern "C" fn mainSimpleSort(
    ptr: *mut u32,
    block: *mut u8,
    quadrant: *mut u16,
    nblock: i32,
    lo: i32,
    hi: i32,
    d: i32,
    budget: *mut i32,
) {
    let ptr = unsafe { from_raw_parts_mut(ptr, (nblock + 8) as usize) };
    let block = unsafe { from_raw_parts_mut(block, (nblock + 8) as usize) };
    let quadrant = unsafe { from_raw_parts_mut(quadrant, (nblock + 8) as usize) };

    main_simple_sort(ptr, block, quadrant, nblock, lo, hi, d, budget);
}

fn main_simple_sort(
    ptr: &mut [u32],
    block: &mut [u8],
    quadrant: &mut [u16],
    nblock: i32,
    lo: i32,
    hi: i32,
    d: i32,
    budget: *mut i32,
) {
    let bigN = hi - lo + 1;
    if bigN < 2 {
        return;
    }

    let mut hp = 0_i32;
    while incs[hp as usize] < bigN {
        hp += 1;
    }

    hp -= 1;

    while hp >= 0 {
        let h = incs[hp as usize];

        let mut i = lo + h;
        loop {
            /*-- copy 1 --*/
            if i > hi {
                break;
            }
            let mut v = ptr[i as usize];
            let mut j = i;
            while main_gtu(
                (ptr[(j - h) as usize] as i32 + d) as u32,
                (v as i32 + d) as u32,
                block,
                quadrant,
                nblock as u32,
                budget,
            ) {
                ptr[j as usize] = ptr[(j - h) as usize];
                j = j - h;
                if j <= (lo + h - 1) {
                    break;
                }
            }
            ptr[j as usize] = v;
            i += 1;

            /*-- copy 2 --*/
            if i > hi {
                break;
            }
            v = ptr[i as usize];
            j = i;
            while main_gtu(
                (ptr[(j - h) as usize] as i32 + d) as u32,
                (v as i32 + d) as u32,
                block,
                quadrant,
                nblock as u32,
                budget,
            ) {
                ptr[j as usize] = ptr[(j - h) as usize];
                j = j - h;
                if j <= (lo + h - 1) {
                    break;
                }
            }
            ptr[j as usize] = v;
            i += 1;

            /*-- copy 3 --*/
            if i > hi {
                break;
            }
            v = ptr[i as usize];
            j = i;
            while main_gtu(
                (ptr[(j - h) as usize] as i32 + d) as u32,
                (v as i32 + d) as u32,
                block,
                quadrant,
                nblock as u32,
                budget,
            ) {
                ptr[j as usize] = ptr[(j - h) as usize];
                j = j - h;
                if j <= (lo + h - 1) {
                    break;
                }
            }
            ptr[j as usize] = v;
            i += 1;

            if unsafe { *budget } < 0 {
                return;
            }
            hp -= 1;
        }
    }
}

#[no_mangle]
pub extern "C" fn BZ2_blockSort(s: &mut EState) {
    let nblock = s.nblock;
    let ptr_arr =
        unsafe { from_raw_parts_mut(s.ptr, (nblock + BZ_N_OVERSHOOT as i32 + 2) as usize) };
    let block_arr =
        unsafe { from_raw_parts_mut(s.block, (nblock + BZ_N_OVERSHOOT as i32 + 2) as usize) };

    let ptr = s.ptr;
    let block = s.block;
    let ftab = s.ftab;

    let verb = s.verbosity;
    let mut wfact = s.workFactor;

    if nblock < 10000 {
        unsafe {
            fallbackSort(s.arr1, s.arr2, ftab, nblock, verb);
        }
    } else {
        // Calculate the location for quadrant, remembering to get
        // the alignment right.  Assumes that &(block[0]) is at least
        // 2-byte aligned -- this should be ok since block is really
        // the first section of arr2.

        let mut i = nblock + BZ_N_OVERSHOOT as i32;
        if (i & 1) != 0 {
            i += 1;
        }
        let quadrant = &mut (block_arr[i as usize]) as *mut u8 as *mut u16;

        // (wfact-1) / 3 puts the default-factor-30
        // transition point at very roughly the same place as
        // with v0.1 and v0.9.0.
        // Not that it particularly matters any more, since the
        // resulting compressed stream is now the same regardless
        // of whether or not we use the main sort or fallback sort.

        if wfact < 1 {
            wfact = 1;
        }

        if wfact > 100 {
            wfact = 100;
        }
        let budgetInit = nblock * ((wfact - 1) / 3);
        let mut budget = budgetInit;

        unsafe {
            mainSort(
                ptr,
                block,
                quadrant,
                ftab,
                nblock,
                verb,
                &mut budget as *mut i32 as *mut u32,
            )
        };

        if verb >= 3 {
            println!(
                "      %{} work, %{} block, ratio %5.2{}\n",
                budgetInit - budget as i32,
                nblock,
                (budgetInit - budget as i32) as f64 / (if nblock == 0 { 1 } else { nblock }) as f64
            );
        }
        if budget < 0 {
            if verb >= 2 {
                println!("    too repetitive; using fallback sorting algorithm");
            }
            unsafe {
                fallbackSort(s.arr1, s.arr2, ftab, nblock, verb);
            }
        }
    }

    s.origPtr = -1;
    for i in 0..s.nblock {
        if ptr_arr[i as usize] == 0 {
            s.origPtr = i;
            break;
        }
    }

    asserth(s.origPtr != -1, 1003);
}

#[no_mangle]
pub extern "C" fn fallbackQSort3(fmap: *mut u32, eclass: *mut u32, loSt: i32, hiSt: i32) {
    let mut r: i32 = 0;
    let mut stackLo: [i32; 100] = [0; 100];
    let mut stackHi: [i32; 100] = [0; 100];

    let fmap = unsafe { from_raw_parts_mut(fmap, (hiSt + 1) as usize) };
    let max_index = fmap
        .iter()
        .skip(loSt as usize)
        .max()
        .cloned()
        .unwrap_or_default() as usize;
    let eclass = unsafe { from_raw_parts_mut(eclass, max_index + 1) };

    stackLo[0] = loSt;
    stackHi[0] = hiSt;
    let mut sp = 1;

    while sp > 0 {
        asserth(sp < 100 - 1, 1004);
        sp -= 1;
        let lo = stackLo[sp as usize];
        let hi = stackHi[sp as usize];

        if hi - lo < 10 {
            fallback_simple_sort(fmap, eclass, lo, hi);
        } else {
            // Random partitioning.  Median of 3 sometimes fails to
            // avoid bad cases.  Median of 9 seems to help but
            // looks rather expensive.  This too seems to work but
            // is cheaper.  Guidance for the magic constants
            // 7621 and 32768 is taken from Sedgewick's algorithms
            // book, chapter 35.
            r = r.wrapping_mul(7621).wrapping_add(1).wrapping_rem(32768);
            let r3 = r.wrapping_rem(3);
            let med = if r3 == 0 {
                eclass[fmap[lo as usize] as usize]
            } else if r3 == 1 {
                eclass[fmap[(lo + hi >> 1) as usize] as usize]
            } else {
                eclass[fmap[hi as usize] as usize]
            };
            let mut ltLo = lo;
            let mut unLo = ltLo;
            let mut gtHi = hi;
            let mut unHi = gtHi;
            loop {
                while !(unLo > unHi) {
                    let n = eclass[fmap[unLo as usize] as usize] as i32 - med as i32;
                    if n == 0 {
                        fmap.swap(unLo as usize, ltLo as usize);
                        ltLo += 1;
                        unLo += 1
                    } else {
                        if n > 0 {
                            break;
                        }
                        unLo += 1
                    }
                }
                while !(unLo > unHi) {
                    let n = eclass[fmap[unHi as usize] as usize] as i32 - med as i32;
                    if n == 0 {
                        fmap.swap(unHi as usize, gtHi as usize);
                        gtHi -= 1;
                        unHi -= 1;
                        continue;
                    } else {
                        if n < 0 {
                            break;
                        }
                        unHi -= 1;
                    }
                }
                if unLo > unHi {
                    break;
                }
                fmap.swap(unLo as usize, unHi as usize);
                unLo += 1;
                unHi -= 1
            }
            assertd(unHi == unLo - 1, "fallbackQSort3(2)\x00");
            if gtHi < ltLo {
                continue;
            }
            let mut n = (ltLo - lo).min(unLo - ltLo);

            fvswap(lo, unLo - n, n, fmap);
            let mut m = (hi - gtHi).min(gtHi - unHi);
            fvswap(unLo, hi - m + 1, m, fmap);
            n = lo + unLo - ltLo - 1;
            m = hi - (gtHi - unHi) + 1;
            if n - lo > hi - m {
                stackLo[sp as usize] = lo;
                stackHi[sp as usize] = n;
                sp += 1;
                stackLo[sp as usize] = m;
                stackHi[sp as usize] = hi;
                sp += 1;
            } else {
                stackLo[sp as usize] = m;
                stackHi[sp as usize] = hi;
                sp += 1;
                stackLo[sp as usize] = lo;
                stackHi[sp as usize] = n;
                sp += 1;
            }
        }
    }
}

fn fvswap(zzp1: i32, zzp2: i32, zzn: i32, fmap: &mut [u32]) {
    let mut yyp1 = zzp1;
    let mut yyp2 = zzp2;
    let mut yyn = zzn;
    while yyn > 0 {
        fmap.swap(yyp1 as usize, yyp2 as usize);
        yyp1 += 1;
        yyp2 += 1;
        yyn -= 1;
    }
}

fn fpush(stackLo: &mut [i32; 100], stackHi: &mut [i32; 100], mut sp: i32, lz: i32, hz: i32) -> i32 {
    stackLo[sp as usize] = lz;
    stackHi[sp as usize] = hz;
    sp += 1;
    sp
}
