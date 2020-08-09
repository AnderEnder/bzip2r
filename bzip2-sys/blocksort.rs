use crate::private_ffi::{mainQSort3, EState, BZ_N_OVERSHOOT, BZ_N_RADIX};
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
    i1: u32,
    i2: u32,
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
        fallbackSort(s.arr1, s.arr2, ftab, nblock, verb);
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

        mainSort(
            ptr,
            block,
            quadrant,
            ftab,
            nblock,
            verb,
            &mut budget as *mut i32,
        );

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
            fallbackSort(s.arr1, s.arr2, ftab, nblock, verb);
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

const FALLBACK_QSORT_SMALL_THRESH: i32 = 10;
const FALLBACK_QSORT_STACK_SIZE: usize = 100;

#[no_mangle]
pub extern "C" fn fallbackQSort3(fmap: *mut u32, eclass: *mut u32, loSt: i32, hiSt: i32) {
    let fmap = unsafe { from_raw_parts_mut(fmap, (hiSt + 1) as usize) };
    let max_index = fmap
        .iter()
        .skip(loSt as usize)
        .max()
        .cloned()
        .unwrap_or_default() as usize;
    let eclass = unsafe { from_raw_parts_mut(eclass, max_index + 1) };
    fallback_qsort3(fmap, eclass, loSt, hiSt)
}

fn fallback_qsort3(fmap: &mut [u32], eclass: &mut [u32], loSt: i32, hiSt: i32) {
    let mut r: i32 = 0;
    let mut stackLo: [i32; FALLBACK_QSORT_STACK_SIZE] = [0; 100];
    let mut stackHi: [i32; FALLBACK_QSORT_STACK_SIZE] = [0; 100];

    stackLo[0] = loSt;
    stackHi[0] = hiSt;
    let mut sp = 1;

    while sp > 0 {
        asserth(sp < FALLBACK_QSORT_STACK_SIZE - 1, 1004);
        sp -= 1;
        let lo = stackLo[sp as usize];
        let hi = stackHi[sp as usize];

        if hi - lo < FALLBACK_QSORT_SMALL_THRESH {
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

#[no_mangle]
pub extern "C" fn fallbackSort(
    fmap: *mut u32,
    eclass: *mut u32,
    bhtab: *mut u32,
    nblock: i32,
    verb: i32,
) {
    let mut ftab: [libc::c_int; 257] = [0; 257];
    let mut ftabCopy: [libc::c_int; 256] = [0; 256];
    let eclass8 = eclass as *mut u8;

    let eclass = unsafe { from_raw_parts_mut(eclass, (nblock + 1) as usize) };
    let eclass8 = unsafe { from_raw_parts_mut(eclass8, (nblock + 1) as usize) };
    let fmap = unsafe { from_raw_parts_mut(fmap, (nblock + 1) as usize) };
    let bhtab = unsafe { from_raw_parts_mut(bhtab, (nblock + 2) as usize) };
    // Initial 1-char radix sort to generate
    // initial fmap and initial BH bits.
    if verb >= 4 {
        println!("        bucket sorting ...\n\x00");
    }

    for i in 0..nblock {
        ftab[eclass8[i as usize] as usize] += 1;
    }
    for i in 0..256 {
        ftabCopy[i as usize] = ftab[i as usize];
    }
    for i in 1..257 {
        ftab[i as usize] += ftab[(i - 1) as usize];
    }
    for i in 0..nblock {
        let j = eclass8[i as usize] as i32;
        let k = ftab[j as usize] - 1;
        ftab[j as usize] = k;
        fmap[k as usize] = i as u32;
    }

    let nBhtab = 2 + nblock / 32;

    for i in 0..nBhtab {
        bhtab[i as usize] = 0;
    }
    for i in 0..256 {
        bhtab[(ftab[i as usize] >> 5 as i32) as usize] |= (1_u32) << (ftab[i as usize] & 31);
    }
    // Inductively refine the buckets.  Kind-of an
    // "exponential radix sort" (!), inspired by the
    // Manber-Myers suffix array construction algorithm.
    // set sentinel bits for block-end detection
    for i in 0..32 {
        bhtab[((nblock + 2 * i) >> 5) as usize] |= (1_u32) << ((nblock + 2 * i) & 31);
        bhtab[((nblock + 2 * i + 1) >> 5) as usize] &= !((1) << ((nblock + 2 * i + 1) & 31));
    }
    // the log(N) loop
    let mut H = 1;
    loop {
        if verb >= 4 {
            println!("        depth %6{} has \x00", H);
        }
        let mut j = 0;
        for i in 0..nblock {
            if bhtab[(i >> 5) as usize] & 1 << (i & 31) != 0 {
                j = i
            }
            let mut k = fmap[i as usize].wrapping_sub(H as u32) as i32;
            if k < 0 {
                k += nblock
            }
            eclass[k as usize] = j as u32;
        }
        let mut nNotDone = 0;
        let mut r = -1;
        loop {
            // find the next non-singleton bucket
            let mut k = r + 1;
            while bhtab[(k >> 5) as usize] & 1_u32 << (k & 31) != 0 && k & 0x1f as i32 != 0 {
                k += 1;
            }
            if bhtab[(k >> 5) as usize] & 1 << (k & 31) != 0 {
                while bhtab[(k >> 5) as usize] == 0xffffffff {
                    k += 32
                }
                while bhtab[(k >> 5) as usize] & (1) << (k & 31) != 0 {
                    k += 1
                }
            }
            let l = k - 1;
            if l >= nblock {
                break;
            }
            while bhtab[(k >> 5) as usize] & 1 << (k & 31) == 0 && k & 0x1f != 0 {
                k += 1
            }
            if bhtab[(k >> 5) as usize] & 1 << (k & 31) == 0 {
                while bhtab[(k >> 5) as usize] == 0 {
                    k += 32
                }
                while bhtab[(k >> 5) as usize] & 1 << (k & 31) == 0 {
                    k += 1
                }
            }
            r = k - 1;
            if r >= nblock {
                break;
            }
            // now [l, r] bracket current bucket
            if r > l {
                nNotDone += r - l + 1;
                fallback_qsort3(fmap, eclass, l, r);
                // scan bucket and generate header bits
                let mut cc = -1;
                for i in l..r + 1 {
                    let cc1 = eclass[fmap[i as usize] as usize] as i32;
                    if cc != cc1 {
                        bhtab[(i >> 5) as usize] |= 1 << (i & 31);
                        cc = cc1
                    }
                }
            }
        }
        if verb >= 4 {
            println!("%6{} unresolved strings\n\x00", nNotDone);
        }
        H *= 2;
        if H > nblock || nNotDone == 0 {
            break;
        }
    }
    // Reconstruct the original block in
    // eclass8 [0 .. nblock-1], since the
    // previous phase destroyed it.
    if verb >= 4 {
        println!("        reconstructing block ...\n\x00");
    }
    let mut j = 0;
    for i in 0..nblock {
        while ftabCopy[j as usize] == 0 {
            j += 1
        }
        ftabCopy[j as usize] -= 1;
        eclass8[fmap[i as usize] as usize] = j as u8;
    }
    asserth(j < 256, 1005);
}

fn bigfreq(b: i32, ftab: &mut [u32]) -> u32 {
    ftab[((b + 1) << 8) as usize] - ftab[(b << 8) as usize]
}

const SETMASK: u32 = 1 << 21;
const CLEARMASK: u32 = !SETMASK;

#[no_mangle]
pub extern "C" fn mainSort(
    ptr_raw: *mut u32,
    block_raw: *mut u8,
    quadrant_raw: *mut u16,
    ftab: *mut u32,
    nblock: i32,
    verb: i32,
    budget: *mut i32,
) {
    let ptr = unsafe { from_raw_parts_mut(ptr_raw, (nblock + BZ_N_OVERSHOOT as i32) as usize) };
    let block = unsafe { from_raw_parts_mut(block_raw, (nblock + BZ_N_OVERSHOOT as i32) as usize) };
    let ftab = unsafe { from_raw_parts_mut(ftab, (nblock + BZ_N_OVERSHOOT as i32) as usize) };
    let quadrant =
        unsafe { from_raw_parts_mut(quadrant_raw, (nblock + BZ_N_OVERSHOOT as i32) as usize) };
    let mut bigDone = [false; 256];
    let mut runningOrder = [0_i32; 256];
    let mut copyStart = [0_i32; 256];
    let mut copyEnd = [0_i32; 256];

    if verb >= 4 {
        println!("        main sort initialise ...");
    }

    // set up the 2-byte frequency table
    for i in 0..65536 + 1 {
        ftab[i] = 0;
    }

    let mut j = (block[0] as i32) << 8;
    let mut i = nblock - 1;
    while i >= 3 {
        quadrant[i as usize] = 0;
        j = (j >> 8) | ((block[i as usize] as i32) << 8);
        ftab[j as usize] += 1;
        quadrant[(i - 1) as usize] = 0;
        j = (j >> 8) | ((block[(i - 1) as usize] as i32) << 8);
        ftab[j as usize] += 1;
        quadrant[(i - 2) as usize] = 0;
        j = (j >> 8) | ((block[(i - 2) as usize] as i32) << 8);
        ftab[j as usize] += 1;
        quadrant[(i - 3) as usize] = 0;
        j = (j >> 8) | ((block[(i - 3) as usize] as i32) << 8);
        ftab[j as usize] += 1;
        i -= 4
    }
    while i >= 0 {
        quadrant[i as usize] = 0;
        j = (j >> 8) | ((block[i as usize] as i32) << 8);
        ftab[j as usize] += 1;
        i -= 1
    }

    // (emphasises close relationship of block & quadrant)
    for i in 0..BZ_N_OVERSHOOT as i32 {
        block[(nblock + i) as usize] = block[i as usize];
        quadrant[(nblock + i) as usize] = 0;
    }

    if verb >= 4 {
        println!("        bucket sorting ...");
    }

    // Complete the initial radix sort
    for i in 1..65536 + 1 {
        ftab[i] += ftab[i - 1];
    }

    let mut s = (block[0] as u16) << 8;
    i = nblock - 1;
    while i >= 3 {
        s = (s >> 8) as u16 | ((block[i as usize] as u16) << 8);
        let mut j = ftab[s as usize] - 1;
        ftab[s as usize] = j;
        ptr[j as usize] = i as u32;
        s = (s >> 8) | ((block[(i - 1) as usize] as u16) << 8);
        j = ftab[s as usize] - 1;
        ftab[s as usize] = j;
        ptr[j as usize] = (i - 1) as u32;
        s = (s >> 8) | ((block[(i - 2) as usize] as u16) << 8);
        j = ftab[s as usize] - 1;
        ftab[s as usize] = j;
        ptr[j as usize] = (i - 2) as u32;
        s = (s >> 8) | ((block[(i - 3) as usize] as u16) << 8);
        j = ftab[s as usize] - 1;
        ftab[s as usize] = j;
        ptr[j as usize] = (i - 3) as u32;
        i -= 4
    }
    while i >= 0 {
        s = (s >> 8) | ((block[i as usize] as u16) << 8);
        let j = ftab[s as usize] - 1;
        ftab[s as usize] = j;
        ptr[j as usize] = i as u32;
        i -= 1;
    }

    // Now ftab contains the first loc of every small bucket.
    // Calculate the running order, from smallest to largest
    // big bucket.

    for i in 0..256 {
        // should be false
        // bigDone[i as usize] = false;
        runningOrder[i as usize] = i;
    }

    {
        let mut h = 1_i32;
        while h <= 256 {
            h = 3 * h + 1;
        }

        loop {
            h = h / 3;
            for i in h..256 {
                let vv = runningOrder[i as usize];
                let mut j = i;
                while bigfreq(runningOrder[(j - h) as usize], ftab) > bigfreq(vv, ftab) {
                    runningOrder[j as usize] = runningOrder[(j - h) as usize];
                    j -= h;
                    if j <= h - 1 {
                        break;
                    }
                }
                runningOrder[j as usize] = vv;
            }
            if h == 1 {
                break;
            }
        }
    }

    // The main sorting loop.
    let mut numQSorted = 0;

    for i in 0..256 {
        //  Process big buckets, starting with the least full.
        //  Basically this is a 3-step process in which we call
        //  mainQSort3 to sort the small buckets [ss, j], but
        //  also make a big effort to avoid the calls if we can.
        let ss = runningOrder[i];

        //  Step 1:
        //  Complete the big bucket [ss] by quicksorting
        //  any unsorted small buckets [ss, j], for j != ss.
        //  Hopefully previous pointer-scanning phases have already
        //  completed many of the small buckets [ss, j], so
        //  we don't have to sort them at all.

        for j in 0..256 {
            if j != ss {
                let sb = (ss << 8) + j;
                if !(ftab[sb as usize] & SETMASK) > 0 {
                    let lo = (ftab[sb as usize] & CLEARMASK) as i32;
                    let hi = (ftab[(sb + 1) as usize] & CLEARMASK) as i32 - 1;
                    if hi > lo {
                        if verb >= 4 {
                            println!(
                                "        qsort [0x%x{}, 0x%x{}]   done %{}   this %{}",
                                ss,
                                j,
                                numQSorted,
                                hi - lo + 1
                            );
                        }
                        unsafe {
                            mainQSort3(
                                ptr_raw,
                                block_raw,
                                quadrant_raw,
                                nblock,
                                lo,
                                hi,
                                BZ_N_RADIX as i32,
                                budget as *mut u32,
                            );
                        }
                        numQSorted += hi - lo + 1;
                        if unsafe { *budget } < 0 {
                            return;
                        }
                    }
                }
                ftab[sb as usize] |= SETMASK;
            }
        }

        asserth(!bigDone[ss as usize], 1006);

        //  Step 2:
        //  Now scan this big bucket [ss] so as to synthesise the
        //  sorted order for small buckets [t, ss] for all t,
        //  including, magically, the bucket [ss,ss] too.
        //  This will avoid doing Real Work in subsequent Step 1's.
        {
            for j in 0..256 {
                copyStart[j as usize] = (ftab[((j << 8) + ss) as usize] & CLEARMASK) as i32;
                copyEnd[j as usize] = (ftab[((j << 8) + ss + 1) as usize] & CLEARMASK) as i32 - 1;
            }

            let start1 = ftab[(ss << 8) as usize] & CLEARMASK;
            let end1 = copyStart[ss as usize];

            for j in start1 as i32..end1 {
                let mut k = ptr[j as usize] as i32 - 1;
                if k < 0 {
                    k += nblock;
                }
                let c1 = block[k as usize];
                if !bigDone[c1 as usize] {
                    copyStart[c1 as usize] += 1;
                    ptr[copyStart[c1 as usize] as usize] = k as u32;
                }
            }

            let start2 = (copyEnd[ss as usize] + 1) as usize;
            let end2 = (ftab[((ss + 1) << 8) as usize] & CLEARMASK) as usize;

            for j in (start2..end2).rev() {
                let mut k = ptr[j as usize] as i32 - 1;
                if k < 0 {
                    k += nblock;
                }
                let c1 = block[k as usize] as usize;
                if !bigDone[c1] {
                    copyEnd[c1] -= 1;
                    ptr[copyEnd[c1] as usize] = k as u32;
                }
            }
        }

        asserth(
            (copyStart[ss as usize] - 1 == copyEnd[ss as usize]) ||
                // Extremely rare case missing in bzip2-1.0.0 and 1.0.1.
                // Necessity for this case is demonstrated by compressing
                // a sequence of approximately 48.5 million of character
                // 251; 1.0.0/1.0.1 will then die here.
                  (copyStart[ss as usize] == 0 && copyEnd[ss as usize] == nblock - 1),
            1007,
        );

        for j in 0..256 {
            ftab[((j << 8) + ss) as usize] |= SETMASK;
        }

        //  Step 3:
        //  The [ss] big bucket is now done.  Record this fact,
        //  and update the quadrant descriptors.  Remember to
        //  update quadrants in the overshoot area too, if
        //  necessary.  The "if (i < 255)" test merely skips
        //  this updating for the last bucket processed, since
        //  updating for the last bucket is pointless.

        //  The quadrant array provides a way to incrementally
        //  cache sort orderings, as they appear, so as to
        //  make subsequent comparisons in fullGtU() complete
        //  faster.  For repetitive blocks this makes a big
        //  difference (but not big enough to be able to avoid
        //  the fallback sorting mechanism, exponential radix sort).

        //  The precise meaning is: at all times:

        //     for 0 <= i < nblock and 0 <= j <= nblock

        //     if block[i] != block[j],

        //        then the relative values of quadrant[i] and
        //             quadrant[j] are meaningless.

        //        else {
        //           if quadrant[i] < quadrant[j]
        //              then the string starting at i lexicographically
        //              precedes the string starting at j

        //           else if quadrant[i] > quadrant[j]
        //              then the string starting at j lexicographically
        //              precedes the string starting at i

        //           else
        //              the relative ordering of the strings starting
        //              at i and j has not yet been determined.
        //        }

        bigDone[ss as usize] = true;

        if i < 255 {
            let bbStart = (ftab[(ss << 8) as usize] & CLEARMASK) as i32;
            let bbSize = (ftab[((ss + 1) << 8) as usize] & CLEARMASK) as i32 - bbStart;
            let mut shifts = 0_i32;

            while (bbSize >> shifts) > 65534 {
                shifts += 1;
            }

            for j in (0..bbSize).rev() {
                let a2update = ptr[(bbStart + j as i32) as usize];
                let qVal = (j >> shifts) as u16;
                quadrant[a2update as usize] = qVal;
                if a2update < BZ_N_OVERSHOOT {
                    quadrant[(a2update as i32 + nblock) as usize] = qVal;
                }
            }
            asserth(((bbSize - 1) >> shifts) <= 65535, 1002);
        }
    }

    if verb >= 4 {
        println!(
            "        %{} pointers, %{} sorted, %{} scanned",
            nblock,
            numQSorted,
            nblock - numQSorted
        );
    }
}
