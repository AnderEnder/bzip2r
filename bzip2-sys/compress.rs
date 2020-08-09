use crate::huffman::BZ2_hbMakeCodeLengths;
use crate::private_ffi::{
    sendMTFValues, BZ2_blockSort, BZ_HDR_h, EState, BZ_G_SIZE, BZ_HDR_0, BZ_HDR_B, BZ_HDR_Z,
    BZ_MAX_SELECTORS, BZ_N_GROUPS, BZ_N_ITERS, BZ_RUNA, BZ_RUNB,
};
use std::slice::from_raw_parts_mut;

#[no_mangle]
pub extern "C" fn BZ2_bsInitWrite(s: &mut EState) {
    s.bsLive = 0;
    s.bsBuff = 0;
}

#[no_mangle]
pub extern "C" fn bsFinishWrite(s: &mut EState) {
    let zbits = unsafe { from_raw_parts_mut(s.zbits, s.nblockMAX as usize) };

    while s.bsLive > 0 {
        zbits[s.numZ as usize] = (s.bsBuff >> 24) as u8;
        s.numZ += 1;
        s.bsBuff <<= 8;
        s.bsLive -= 8;
    }
}

fn bsNEEDW(s: &mut EState) {
    let zbits = unsafe { from_raw_parts_mut(s.zbits, s.nblockMAX as usize) };

    while s.bsLive >= 8 {
        zbits[s.numZ as usize] = (s.bsBuff >> 24) as u8;
        s.numZ += 1;
        s.bsBuff <<= 8;
        s.bsLive -= 8;
    }
}

#[no_mangle]
pub extern "C" fn bsW(s: &mut EState, n: i32, v: u32) {
    bsNEEDW(s);
    s.bsBuff |= v << (32 - s.bsLive - n);
    s.bsLive += n;
}

#[no_mangle]
pub extern "C" fn bsPutUInt32(s: &mut EState, u: u32) {
    bsW(s, 8, (u >> 24) & 0xff);
    bsW(s, 8, (u >> 16) & 0xff);
    bsW(s, 8, (u >> 8) & 0xff);
    bsW(s, 8, u & 0xff);
}

#[no_mangle]
pub extern "C" fn bsPutUChar(s: &mut EState, c: u8) {
    bsW(s, 8, c as u32);
}

#[no_mangle]
pub extern "C" fn makeMaps_e(s: &mut EState) {
    s.nInUse = 0;

    for (in_use, unseq_to_seq) in s.inUse.iter().zip(s.unseqToSeq.iter_mut()) {
        if *in_use == 1 {
            *unseq_to_seq = s.nInUse as u8;
            s.nInUse += 1;
        }
    }
}

pub fn assertd(_cond: bool, _msg: &'static str) {
    1;
}

pub fn asserth(_cond: bool, _msg: i32) {
    1;
}

#[no_mangle]
pub extern "C" fn generateMTFValues(s: &mut EState) {
    /*
       After sorting (eg, here),
          s->arr1 [ 0 .. s->nblock-1 ] holds sorted order,
          and
          ((UChar*)s->arr2) [ 0 .. s->nblock-1 ]
          holds the original block data.

       The first thing to do is generate the MTF values,
       and put them in
          ((UInt16*)s->arr1) [ 0 .. s->nblock-1 ].
       Because there are strictly fewer or equal MTF values
       than block values, ptr values in this area are overwritten
       with MTF values only when they are no longer needed.

       The final compressed bitstream is generated into the
       area starting at
          (UChar*) (&((UChar*)s->arr2)[s->nblock])

       These storage aliases are set up in bzCompressInit(),
       except for the last one, which is arranged in
       compressBlock().
    */
    let ptr = unsafe { from_raw_parts_mut(s.ptr, (s.nblockMAX + 2) as usize) };
    let block = unsafe { from_raw_parts_mut(s.block, (s.nblockMAX + 2) as usize) };
    let mtfv = unsafe { from_raw_parts_mut(s.mtfv, (s.nblockMAX + 2) as usize) };

    makeMaps_e(s);
    let EOB = s.nInUse + 1;

    for i in 0..(EOB + 1) as usize {
        s.mtfFreq[i] = 0;
    }

    let mut wr = 0;
    let mut zPend = 0;
    let mut yy = [0_u8; 256];

    for i in 0..s.nInUse as usize {
        yy[i] = i as u8;
    }

    for i in 0..s.nblock as usize {
        // assertd(wr <= i, "generateMTFValues(1)");

        let mut j = ptr[i] as i32 - 1;
        if j < 0 {
            j += s.nblock;
        }

        let ll_i = s.unseqToSeq[block[j as usize] as usize];
        // assertd((ll_i as i32) < s.nInUse, "generateMTFValues(2a)");

        if yy[0] == ll_i {
            zPend += 1;
        } else {
            if zPend > 0 {
                zPend -= 1;
                loop {
                    if (zPend & 1) != 0 {
                        mtfv[wr] = BZ_RUNB as u16;
                        wr += 1;
                        s.mtfFreq[BZ_RUNB as usize] += 1;
                    } else {
                        mtfv[wr] = BZ_RUNA as u16;
                        wr += 1;
                        s.mtfFreq[BZ_RUNA as usize] += 1;
                    }
                    if zPend < 2 {
                        break;
                    }
                    zPend = (zPend - 2) / 2;
                }
                zPend = 0;
            }
            {
                let mut rtmp = yy[1];
                yy[1] = yy[0];
                let rll_i = ll_i;
                let mut ryy_j = 1;

                while rll_i != rtmp {
                    ryy_j += 1;
                    let rtmp2 = rtmp;
                    rtmp = yy[ryy_j];
                    yy[ryy_j] = rtmp2;
                }

                yy[0] = rtmp;
                //j = yy[1] as i32 - yy[0] as i32;
                j = ryy_j as i32;
                mtfv[wr] = (j + 1) as u16;
                wr += 1;
                s.mtfFreq[(j + 1) as usize] += 1;
            }
        }
    }

    if zPend > 0 {
        zPend -= 1;
        loop {
            if (zPend & 1) != 0 {
                mtfv[wr] = BZ_RUNB as u16;
                wr += 1;
                s.mtfFreq[BZ_RUNB as usize] += 1;
            } else {
                mtfv[wr] = BZ_RUNA as u16;
                wr += 1;
                s.mtfFreq[BZ_RUNA as usize] += 1;
            }
            if zPend < 2 {
                break;
            }
            zPend = (zPend - 2) / 2;
        }
        //zPend = 0;
    }

    mtfv[wr] = EOB as u16;
    wr += 1;
    s.mtfFreq[EOB as usize] += 1;
    s.nMTF = wr as i32;
}

fn BZ_ITER_50(s: &mut EState, mtfv: &mut [u16], gs: usize) -> (u32, u32, u32) {
    let mut cost01 = 0;
    let mut cost23 = 0;
    let mut cost45 = 0;

    for nn in 0..50 {
        let icv = mtfv[gs + nn];
        cost01 += s.len_pack[icv as usize][0] as u32;
        cost23 += s.len_pack[icv as usize][1] as u32;
        cost45 += s.len_pack[icv as usize][2] as u32;
    }
    (cost01, cost23, cost45)
}

fn BZ_ITUR_50(s: &mut EState, mtfv: &mut [u16], bt: usize, gs: usize) {
    for nn in 0..50 {
        s.rfreq[bt][mtfv[gs + nn] as usize] += 1
    }
}

const BZ_GREATER_ICOST: u8 = 15;
const BZ_LESSER_ICOST: u8 = 0;

#[no_mangle]
pub extern "C" fn sendMTFValues2(s: &mut EState) {
    //    Int32 v, t, i, j, gs, ge, totc, bt, bc, iter;
    //    Int32 nSelectors, alphaSize, minLen, maxLen, selCtr;
    //    Int32 nGroups, nBytes;
    let mut nSelectors = 0;

    /*--
    UChar  len [BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
    is a global since the decoder also needs it.

    Int32  code[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
    Int32  rfreq[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
    are also globals only used in this proc.
    Made global to keep stack frame size small.
    --*/

    let mut cost = [0_u16; BZ_N_GROUPS as usize];
    let mut fave = [0_i32; BZ_N_GROUPS as usize];
    // let mut gs = 0;
    // let mut ge: i32 = 0;

    let mtfv = unsafe { from_raw_parts_mut(s.mtfv, (s.nblockMAX + 2) as usize) };

    if s.verbosity >= 3 {
        println!(
            "      {} in block, {} after MTF & 1-2 coding,
%{}+2 syms in use",
            s.nblock, s.nMTF, s.nInUse
        );
    }
    let alphaSize = s.nInUse + 2;
    for t in 0..BZ_N_GROUPS as usize {
        for v in 0..alphaSize as usize {
            s.len[t][v] = BZ_GREATER_ICOST;
        }
    }
    /*--- Decide how many coding tables to use ---*/
    asserth(s.nMTF > 0, 3001);

    let nGroups = match s.nMTF {
        m if m < 200 => 2,
        m if m < 600 => 3,
        m if m < 1200 => 4,
        m if m < 2400 => 5,
        _ => 6,
    };

    //*--- Generate an initial set of coding tables ---*/
    {
        // Int32 nPart, remF, tFreq, aFreq;

        let mut nPart = nGroups;
        let mut remF = s.nMTF;
        let mut gs = 0_i32;

        while nPart > 0 {
            let tFreq = remF / nPart as i32;
            let mut ge = gs - 1;
            let mut aFreq = 0;
            while aFreq < tFreq && (ge as i32) < (alphaSize - 1) {
                ge += 1;
                aFreq += s.mtfFreq[ge as usize];
            }

            if ge > gs && nPart != nGroups && nPart != 1 && ((nGroups - nPart) % 2 == 1) {
                aFreq -= s.mtfFreq[ge as usize];
                ge -= 1;
            }

            if s.verbosity >= 3 {
                println!(
                    "      initial group {}, [{} .. {}] has {} syms ({}%%)",
                    nPart,
                    gs,
                    ge,
                    aFreq,
                    (100.0 * aFreq as f64) / s.nMTF as f64
                );
            }

            for v in 0..alphaSize {
                if v >= gs && v <= ge {
                    s.len[(nPart - 1)][v as usize] = BZ_LESSER_ICOST;
                } else {
                    s.len[(nPart - 1)][v as usize] = BZ_GREATER_ICOST;
                }
            }

            nPart -= 1;
            gs = ge + 1;
            remF -= aFreq;
        }
    }

    // Iterate up to BZ_N_ITERS times to improve the tables.
    for iter in 0..BZ_N_ITERS as usize {
        for t in 0..nGroups as usize {
            fave[t] = 0;
        }

        for t in 0..nGroups as usize {
            for v in 0..alphaSize as usize {
                s.rfreq[t][v] = 0;
            }
        }

        // Set up an auxiliary length table which is used to fast-track
        // the common case (nGroups == 6).
        if nGroups == 6 {
            for v in 0..alphaSize as usize {
                s.len_pack[v][0] = ((s.len[1][v] as u32) << 16) | (s.len[0][v] as u32);
                s.len_pack[v][1] = ((s.len[3][v] as u32) << 16) | (s.len[2][v] as u32);
                s.len_pack[v][2] = ((s.len[5][v] as u32) << 16) | (s.len[4][v] as u32);
            }
        }

        nSelectors = 0;
        let mut totc = 0;
        let mut gs = 0;
        loop {
            // Set group start & end marks. --*/
            if gs >= s.nMTF {
                break;
            }

            let mut ge = gs + BZ_G_SIZE as i32 - 1;

            if ge >= s.nMTF {
                ge = s.nMTF - 1;
            }

            // Calculate the cost of this group as coded
            // by each of the coding tables.
            for t in 0..nGroups as usize {
                cost[t] = 0;
            }

            if nGroups == 6 && 50 == ge - gs + 1 {
                //    register UInt32 cost01, cost23, cost45;
                //    register UInt16 icv;
                let (cost01, cost23, cost45) = BZ_ITER_50(s, mtfv, gs as usize);
                cost[0] = cost01 as u16 & 0xffff;
                cost[1] = (cost01 >> 16) as u16;
                cost[2] = cost23 as u16 & 0xffff;
                cost[3] = (cost23 >> 16) as u16;
                cost[4] = cost45 as u16 & 0xffff;
                cost[5] = (cost45 >> 16) as u16;
            } else {
                // slow version which correctly handles all situations ---*/
                for i in gs..ge + 1 {
                    let icv = mtfv[i as usize] as usize;
                    for t in 0..nGroups as usize {
                        cost[t] += s.len[t][icv] as u16;
                    }
                }
            }

            // Find the coding table which is best for this group,
            // and record its identity in the selector table.

            let mut bc = 999999999_i32;
            let mut bt = -1;
            for t in 0..nGroups as usize {
                if (cost[t] as i32) < bc {
                    bc = cost[t] as i32;
                    bt = t as i32;
                };
            }
            totc += bc;
            fave[bt as usize] += 1;
            s.selector[nSelectors] = bt as u8;
            nSelectors += 1;

            //  Increment the symbol frequencies for the selected table.
            if nGroups == 6 && 50 == ge - gs + 1 {
                // fast track the common case ---*/
                // BZ_ITUR_50(s, mtfv, bt as usize, gs as usize);
                for nn in 0..50 {
                    s.rfreq[bt as usize][mtfv[gs as usize + nn] as usize] += 1
                }
            } else {
                // slow version which correctly handles all situations
                for i in gs..ge + 1 {
                    s.rfreq[bt as usize][mtfv[i as usize] as usize] += 1;
                }
            }

            gs = ge + 1;
        }

        if s.verbosity >= 3 {
            println!(
                "      pass {}: size is {}, grp uses are ",
                iter + 1,
                totc / 8,
            );
            for t in 0..nGroups as usize {
                println!("{} ", fave[t]);
                println!("");
            }
        }

        // Recompute the tables based on the accumulated frequencies.
        // maxLen was changed from 20 to 17 in bzip2-1.0.3.  See
        // comment in huffman.c for details
        for t in 0..nGroups as usize {
            BZ2_hbMakeCodeLengths(
                &mut (s.len[t][0]),
                &mut (s.rfreq[t][0]),
                alphaSize,
                17, // 20
            );
        }
    }

    asserth(nGroups < 8, 3002);
    asserth(
        nSelectors < 32768 && nSelectors <= BZ_MAX_SELECTORS as usize,
        3003,
    );

    {
        // UChar pos[BZ_N_GROUPS], ll_i, tmp2, tmp;
        let mut pos = [0; BZ_N_GROUPS as usize];
        for i in 0..nGroups as usize {
            pos[i] = i;
        }

        for i in 0..nSelectors as usize {
            let ll_i = s.selector[i] as usize;
            let mut j = 0;
            let mut tmp = pos[j];
            while ll_i != tmp {
                j += 1;
                let tmp2 = tmp;
                tmp = pos[j];
                pos[j] = tmp2;
            }
            pos[0] = tmp;
            s.selectorMtf[i] = j as u8;
        }
    }
}

#[no_mangle]
pub extern "C" fn BZ2_compressBlock(s: &mut EState, is_last_block: u8) {
    if s.nblock > 0 {
        s.blockCRC = !s.blockCRC;
        s.combinedCRC = (s.combinedCRC << 1) | (s.combinedCRC >> 31);
        s.combinedCRC ^= s.blockCRC;

        if s.blockNo > 1 {
            s.numZ = 0;
        }

        if s.verbosity >= 2 {
            // fixup later
            println!(
                "block %{}: crc = 0x{}, combined CRC = 0x{}, size = %{}",
                s.blockNo, s.blockCRC, s.combinedCRC, s.nblock
            );
        }

        unsafe {
            BZ2_blockSort(s);
        }
    }

    // s.zbits = (UChar *)(&((UChar *)s.arr2)[s.nblock as usize]);
    // unsafe {
    //     let array_raw = s.arr2 as *mut u8;
    //     let bits = from_raw_parts_mut(array_raw, (s.nblock + 1) as usize);
    //     let raw = bits[s.nblock as usize] as *mut u8;
    //     s.zbits = raw;
    // }
    set_zbits(s);

    // /*-- If this is the first block, create the stream header. --*/
    if s.blockNo == 1 {
        BZ2_bsInitWrite(s);
        bsPutUChar(s, BZ_HDR_B as u8);
        bsPutUChar(s, BZ_HDR_Z as u8);
        bsPutUChar(s, BZ_HDR_h as u8);
        bsPutUChar(s, (BZ_HDR_0 as i32 + s.blockSize100k) as u8);
    }

    if s.nblock > 0 {
        bsPutUChar(s, 0x31);
        bsPutUChar(s, 0x41);
        bsPutUChar(s, 0x59);
        bsPutUChar(s, 0x26);
        bsPutUChar(s, 0x53);
        bsPutUChar(s, 0x59);
        // *-- Now the block's CRC, so it is in a known place. --*/
        bsPutUInt32(s, s.blockCRC);
        //*--
        //   Now a single bit indicating (non-)randomisation.
        //   As of version 0.9.5, we use a better sorting algorithm
        //   which makes randomisation unnecessary.  So always set
        //   the randomised bit to 'no'.  Of course, the decoder
        //   still needs to be able to handle randomised blocks
        //   so as to maintain backwards compatibility with
        //   older versions of bzip2.
        //--*/
        bsW(s, 1, 0);

        bsW(s, 24, s.origPtr as u32);
        generateMTFValues(s);
        unsafe {
            sendMTFValues(s);
        }
    }

    //*-- If this is the last block, add the stream trailer. --*/
    if is_last_block == 1 {
        bsPutUChar(s, 0x17);
        bsPutUChar(s, 0x72);
        bsPutUChar(s, 0x45);
        bsPutUChar(s, 0x38);
        bsPutUChar(s, 0x50);
        bsPutUChar(s, 0x90);
        bsPutUInt32(s, s.combinedCRC);
        if s.verbosity >= 2 {
            print!("    final combined CRC = 0x{}\n   ", s.combinedCRC);
        }
        bsFinishWrite(s);
    }
}

fn set_zbits(s: &mut EState) {
    s.zbits = unsafe { &mut *(s.arr2 as *mut u8).offset(s.nblock as isize) as *mut u8 };
}
