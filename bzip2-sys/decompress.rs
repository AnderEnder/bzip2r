use crate::bzlib::BZ2_indexIntoF;
use crate::bzlib::BZALLOC;
use crate::huffman::bz2_hb_create_decode_tables;
use crate::private_ffi::DState;
use crate::private_ffi::{
    BZ_HDR_h, BZ_DATA_ERROR, BZ_DATA_ERROR_MAGIC, BZ_HDR_0, BZ_HDR_B, BZ_HDR_Z, BZ_MEM_ERROR,
    BZ_N_GROUPS, BZ_OK, BZ_X_MAGIC_1, BZ_X_MAGIC_2, BZ_X_MAGIC_3, BZ_X_MAGIC_4, BZ_X_RANDBIT,
};
use crate::private_ffi::{BZ_X_BCRC_1, BZ_X_BCRC_2, BZ_X_BCRC_3, BZ_X_BCRC_4};
use crate::private_ffi::{
    BZ_X_BLKHDR_1, BZ_X_BLKHDR_2, BZ_X_BLKHDR_3, BZ_X_BLKHDR_4, BZ_X_BLKHDR_5, BZ_X_BLKHDR_6,
};
use crate::private_ffi::{BZ_X_ORIGPTR_1, BZ_X_ORIGPTR_2, BZ_X_ORIGPTR_3};
use crate::private_ffi::{MTFA_SIZE, MTFL_SIZE};
use crate::randtable::BZ2_rNums;

use std::mem::size_of;
use std::slice::from_raw_parts_mut;

#[no_mangle]
pub extern "C" fn makeMaps_d(s: &mut DState) {
    s.nInUse = 0;
    for i in 0..256 {
        if s.inUse[i] > 0 {
            s.seqToUnseq[s.nInUse as usize] = i as u8;
            s.nInUse += 1;
        }
    }
}

fn get_bits(s: &mut DState, nnn: i32) -> i32 {
    // s.state = lll;
    let strm = unsafe { s.strm.as_mut() }.unwrap();
    let next_in = unsafe { strm.next_in.as_mut() }.unwrap();
    loop {
        if s.bsLive >= nnn {
            let v = (s.bsBuff >> (s.bsLive - nnn)) & ((1 << nnn) - 1);
            s.bsLive -= nnn;
            return v as i32;
        }
        if strm.avail_in == 0 {
            return BZ_OK as i32;
        }
        s.bsBuff = (s.bsBuff << 8) | (*next_in as u32);
        s.bsLive += 8;

        *next_in += 1;
        strm.avail_in -= 1;
        strm.total_in_lo32 += 1;
        if strm.total_in_lo32 == 0 {
            strm.total_in_hi32 += 1;
        }
    }
}

fn BZ2_decompress(s: &mut DState) -> i32 {
    if s.state == BZ_X_MAGIC_1 as i32 {
        // initialise the save area
        s.save_i = 0;
        s.save_j = 0;
        s.save_t = 0;
        s.save_alphaSize = 0;
        s.save_nGroups = 0;
        s.save_nSelectors = 0;
        s.save_EOB = 0;
        s.save_groupNo = 0;
        s.save_groupPos = 0;
        s.save_nextSym = 0;
        s.save_nblockMAX = 0;
        s.save_nblock = 0;
        s.save_es = 0;
        s.save_N = 0;
        s.save_curr = 0;
        s.save_zt = 0;
        s.save_zn = 0;
        s.save_zvec = 0;
        s.save_zj = 0;
        s.save_gSel = 0;
        s.save_gMinlen = 0;
        s.save_gLimit = std::ptr::null_mut();
        s.save_gBase = std::ptr::null_mut();
        s.save_gPerm = std::ptr::null_mut();
    }

    // restore from the save area
    let i = s.save_i;
    let j = s.save_j;
    let t = s.save_t;
    let alphaSize = s.save_alphaSize;
    let nGroups = s.save_nGroups;
    let nSelectors = s.save_nSelectors;
    let EOB = s.save_EOB;
    let groupNo = s.save_groupNo;
    let groupPos = s.save_groupPos;
    let nextSym = s.save_nextSym;
    let nblockMAX = s.save_nblockMAX;
    let nblock = s.save_nblock;
    let es = s.save_es;
    let N = s.save_N;
    let curr = s.save_curr;
    let zt = s.save_zt;
    let zn = s.save_zn;
    let zvec = s.save_zvec;
    let zj = s.save_zj;
    let gSel = s.save_gSel;
    let gMinlen = s.save_gMinlen;
    let gLimit = s.save_gLimit;
    let gBase = s.save_gBase;
    let gPerm = s.save_gPerm;

    let mut retVal = BZ_OK as i32;
    let strm = unsafe { s.strm.as_mut() }.unwrap();

    // let uc = 0;
    'outer: loop {
        let mut fallthrough = false;
        if s.state as u32 == BZ_X_MAGIC_1 {
            // GET_UCHAR(BZ_X_MAGIC_1, uc);
            // GET_BITS(lll, uuu, 8)
            let uc = get_bits(s, 8);

            if uc == 0 {
                retVal = 0;
                break 'outer;
            } else if uc != BZ_HDR_B as i32 {
                retVal = BZ_DATA_ERROR_MAGIC;
                break 'outer;
            }
            fallthrough = true;
        }
        if fallthrough || s.state as u32 == BZ_X_MAGIC_2 {
            let uc = get_bits(s, 8);

            if uc == 0 {
                retVal = 0;
                break 'outer;
            } else if uc != BZ_HDR_Z as i32 {
                retVal = BZ_DATA_ERROR_MAGIC;
                break 'outer;
            }
            fallthrough = true;
            // GET_UCHAR(BZ_X_MAGIC_2, uc);
            // if (uc != BZ_HDR_Z)
            //    RETURN(BZ_DATA_ERROR_MAGIC);
        }
        if fallthrough || s.state as u32 == BZ_X_MAGIC_3 {
            let uc = get_bits(s, 8);

            if uc == 0 {
                retVal = 0;
                break 'outer;
            } else if uc != BZ_HDR_h as i32 {
                retVal = BZ_DATA_ERROR_MAGIC;
                break 'outer;
            }
            fallthrough = true;
        }
        if fallthrough || s.state as u32 == BZ_X_MAGIC_4 {
            let uc = get_bits(s, 8);

            if uc == 0 {
                retVal = 0;
                break 'outer;
            }
            s.blockSize100k = uc;
            if s.blockSize100k < (BZ_HDR_0 as i32 + 1) || s.blockSize100k > (BZ_HDR_0 as i32 + 9) {
                retVal = BZ_DATA_ERROR_MAGIC;
                break 'outer;
            }
            s.blockSize100k -= BZ_HDR_0 as i32;
            if s.smallDecompress > 0 {
                s.ll16 =
                    BZALLOC(strm, s.blockSize100k * 100000 * size_of::<u16>() as i32) as *mut u16;
                s.ll4 = BZALLOC(
                    strm,
                    ((1 + s.blockSize100k * 100000) >> 1) * size_of::<u8>() as i32,
                ) as *mut u8;
                if s.ll16.is_null() || s.ll4.is_null() {
                    retVal = BZ_MEM_ERROR;
                    break 'outer;
                }
            } else {
                s.tt =
                    BZALLOC(strm, s.blockSize100k * 100000 * size_of::<u32>() as i32) as *mut u32;
                if s.tt.is_null() {
                    retVal = BZ_MEM_ERROR;
                    break 'outer;
                }
            }
            fallthrough = true;
        }
        if fallthrough || s.state as u32 == BZ_X_BLKHDR_1 {
            let uc = get_bits(s, 8);

            if uc == 0 {
                retVal = 0;
                break 'outer;
            } else if uc == 0x17 {
                //goto endhdr_2;
            } else if uc != 0x31 {
                retVal = BZ_DATA_ERROR;
                break 'outer;
            }
            fallthrough = true;
        }
        if fallthrough || s.state as u32 == BZ_X_BLKHDR_2 {
            let uc = get_bits(s, 8);
            if uc == 0 {
                retVal = 0;
                break 'outer;
            } else if uc != 0x41 {
                retVal = BZ_DATA_ERROR;
                break 'outer;
            }
            fallthrough = true;
        }
        if fallthrough || s.state as u32 == BZ_X_BLKHDR_3 {
            //    GET_UCHAR(BZ_X_BLKHDR_3, uc);
            //    if (uc != 0x59)
            //       RETURN(BZ_DATA_ERROR);
            fallthrough = true;
        }
        if fallthrough || s.state as u32 == BZ_X_BLKHDR_4 {
            //    GET_UCHAR(BZ_X_BLKHDR_4, uc);
            //    if (uc != 0x26)
            //       RETURN(BZ_DATA_ERROR);
            fallthrough = true;
        }
        // BZ_X_BLKHDR_5 => {
        //     //    GET_UCHAR(BZ_X_BLKHDR_5, uc);
        //     //    if (uc != 0x53)
        //     //       RETURN(BZ_DATA_ERROR);
        // }
        // BZ_X_BLKHDR_6 => {
        //     //    GET_UCHAR(BZ_X_BLKHDR_6, uc);
        //     //    if (uc != 0x59)
        //     //       RETURN(BZ_DATA_ERROR);

        //     //    s.currBlockNo+=1;
        //     //    if (s.verbosity >= 2)
        //     //       VPrintf1("\n    [%d: huff+mtf ", s.currBlockNo);

        //     //    s.storedBlockCRC = 0;
        // }
        // BZ_X_BCRC_1 => {
        //     //    GET_UCHAR(BZ_X_BCRC_1, uc);
        //     //    s.storedBlockCRC = (s.storedBlockCRC << 8) | ((UInt32)uc);
        // }
        // BZ_X_BCRC_2 => {
        //     //    GET_UCHAR(BZ_X_BCRC_2, uc);
        //     //    s.storedBlockCRC = (s.storedBlockCRC << 8) | ((UInt32)uc);
        // }
        // BZ_X_BCRC_3 => {
        //     //    GET_UCHAR(BZ_X_BCRC_3, uc);
        //     //    s.storedBlockCRC = (s.storedBlockCRC << 8) | ((UInt32)uc);
        // }
        // BZ_X_BCRC_4 => {
        //     //    GET_UCHAR(BZ_X_BCRC_4, uc);
        //     //    s.storedBlockCRC = (s.storedBlockCRC << 8) | ((UInt32)uc);
        // }
        // BZ_X_RANDBIT => {
        //     //    GET_BITS(BZ_X_RANDBIT, s.blockRandomised, 1);
        //     //    s.origPtr = 0;
        // }
        // BZ_X_ORIGPTR_1 => {
        //     //    GET_UCHAR(BZ_X_ORIGPTR_1, uc);
        //     //    s.origPtr = (s.origPtr << 8) | ((Int32)uc);
        // }
        // BZ_X_ORIGPTR_2 => {
        //     //    GET_UCHAR(BZ_X_ORIGPTR_2, uc);
        //     //    s.origPtr = (s.origPtr << 8) | ((Int32)uc);
        // }
        // BZ_X_ORIGPTR_3 => {
        //     //    GET_UCHAR(BZ_X_ORIGPTR_3, uc);
        //     //    s.origPtr = (s.origPtr << 8) | ((Int32)uc);
        //     //    if (s.origPtr < 0)
        //     //       RETURN(BZ_DATA_ERROR);
        //     //    if (s.origPtr > 10 + 100000 * s.blockSize100k)
        //     //       RETURN(BZ_DATA_ERROR);
        // }
        //       /*-=1- Receive the mapping table -=1-*/
        //       for (i = 0; i < 16; i+=1)
        //       {
        //          GET_BIT(BZ_X_MAPPING_1, uc);
        //          if (uc == 1)
        //             s.inUse16[i] = True;
        //          else
        //             s.inUse16[i] = False;
        //       }

        //       for (i = 0; i < 256; i+=1)
        //          s.inUse[i] = False;

        //       for (i = 0; i < 16; i+=1)
        //          if (s.inUse16[i])
        //             for (j = 0; j < 16; j+=1)
        //             {
        //                GET_BIT(BZ_X_MAPPING_2, uc);
        //                if (uc == 1)
        //                   s.inUse[i * 16 + j] = True;
        //             }
        //       makeMaps_d(s);
        //       if (s.nInUse == 0)
        //          RETURN(BZ_DATA_ERROR);
        //       alphaSize = s.nInUse + 2;

        //       /*-=1- Now the selectors -=1-*/
        //       GET_BITS(BZ_X_SELECTOR_1, nGroups, 3);
        //       if (nGroups < 2 || nGroups > BZ_N_GROUPS)
        //          RETURN(BZ_DATA_ERROR);
        //       GET_BITS(BZ_X_SELECTOR_2, nSelectors, 15);
        //       if (nSelectors < 1)
        //          RETURN(BZ_DATA_ERROR);
        //       for (i = 0; i < nSelectors; i+=1)
        //       {
        //          j = 0;
        //          while (True)
        //          {
        //             GET_BIT(BZ_X_SELECTOR_3, uc);
        //             if (uc == 0)
        //                break;
        //             j+=1;
        //             if (j >= nGroups)
        //                RETURN(BZ_DATA_ERROR);
        //          }
        //          /* Having more than BZ_MAX_SELECTORS doesn't make much sense
        //             since they will never be used, but some implementations might
        //             "round up" the number of selectors, so just ignore those. */
        //          if (i < BZ_MAX_SELECTORS)
        //             s.selectorMtf[i] = j;
        //       }
        //       if (nSelectors > BZ_MAX_SELECTORS)
        //          nSelectors = BZ_MAX_SELECTORS;

        //       /*-=1- Undo the MTF values for the selectors. -=1-*/
        //       {
        //          UChar pos[BZ_N_GROUPS], tmp, v;
        //          for (v = 0; v < nGroups; v+=1)
        //             pos[v] = v;

        //          for (i = 0; i < nSelectors; i+=1)
        //          {
        //             v = s.selectorMtf[i];
        //             tmp = pos[v];
        //             while (v > 0)
        //             {
        //                pos[v] = pos[v - 1];
        //                v-=1;
        //             }
        //             pos[0] = tmp;
        //             s.selector[i] = tmp;
        //          }
        //       }

        //       /*-=1- Now the coding tables -=1-*/
        //       for (t = 0; t < nGroups; t+=1)
        //       {
        //          GET_BITS(BZ_X_CODING_1, curr, 5);
        //          for (i = 0; i < alphaSize; i+=1)
        //          {
        //             while (True)
        //             {
        //                if (curr < 1 || curr > 20)
        //                   RETURN(BZ_DATA_ERROR);
        //                GET_BIT(BZ_X_CODING_2, uc);
        //                if (uc == 0)
        //                   break;
        //                GET_BIT(BZ_X_CODING_3, uc);
        //                if (uc == 0)
        //                   curr+=1;
        //                else
        //                   curr-=1;
        //             }
        //             s.len[t][i] = curr;
        //          }
        //       }

        //       /*-=1- Create the Huffman decoding tables -=1-*/
        //       for (t = 0; t < nGroups; t+=1)
        //       {
        //          minLen = 32;
        //          maxLen = 0;
        //          for (i = 0; i < alphaSize; i+=1)
        //          {
        //             if (s.len[t][i] > maxLen)
        //                maxLen = s.len[t][i];
        //             if (s.len[t][i] < minLen)
        //                minLen = s.len[t][i];
        //          }
        //          BZ2_hbCreateDecodeTables(
        //              &(s.limit[t][0]),
        //              &(s.base[t][0]),
        //              &(s.perm[t][0]),
        //              &(s.len[t][0]),
        //              minLen, maxLen, alphaSize);
        //          s.minLens[t] = minLen;
        //       }

        //       /*-=1- Now the MTF values -=1-*/
        //       EOB = s.nInUse + 1;
        //       nblockMAX = 100000 * s.blockSize100k;
        //       groupNo = -1;
        //       groupPos = 0;

        //       for (i = 0; i <= 255; i+=1)
        //          s.unzftab[i] = 0;

        //       /*-=1 MTF init -=1*/
        //       {
        //          Int32 ii, jj, kk;
        //          kk = MTFA_SIZE - 1;
        //          for (ii = 256 / MTFL_SIZE - 1; ii >= 0; ii-=1)
        //          {
        //             for (jj = MTFL_SIZE - 1; jj >= 0; jj-=1)
        //             {
        //                s.mtfa[kk] = (UChar)(ii * MTFL_SIZE + jj);
        //                kk-=1;
        //             }
        //             s.mtfbase[ii] = kk + 1;
        //          }
        //       }
        //       /*-=1 end MTF init -=1*/
        //       nblock = 0;
        //       GET_MTF_VAL(BZ_X_MTF_1, BZ_X_MTF_2, nextSym);

        //       while (True)
        //       {

        //          if (nextSym == EOB)
        //             break;

        //          if (nextSym == BZ_RUNA || nextSym == BZ_RUNB)
        //          {

        //             es = -1;
        //             N = 1;
        //             do
        //             {
        //                /* Check that N doesn't get too big, so that es doesn't
        //                   go negative.  The maximum value that can be
        //                   RUNA/RUNB encoded is equal to the block size (post
        //                   the initial RLE), viz, 900k, so bounding N at 2
        //                   million should guard against overflow without
        //                   rejecting any legitimate inputs. */
        //                if (N >= 2 * 1024 * 1024)
        //                   RETURN(BZ_DATA_ERROR);
        //                if (nextSym == BZ_RUNA)
        //                   es = es + (0 + 1) * N;
        //                else if (nextSym == BZ_RUNB)
        //                   es = es + (1 + 1) * N;
        //                N = N * 2;
        //                GET_MTF_VAL(BZ_X_MTF_3, BZ_X_MTF_4, nextSym);
        //             } while (nextSym == BZ_RUNA || nextSym == BZ_RUNB);

        //             es+=1;
        //             uc = s.seqToUnseq[s.mtfa[s.mtfbase[0]]];
        //             s.unzftab[uc] += es;

        //             if (s.smallDecompress)
        //                while (es > 0)
        //                {
        //                   if (nblock >= nblockMAX)
        //                      RETURN(BZ_DATA_ERROR);
        //                   s.ll16[nblock] = (UInt16)uc;
        //                   nblock+=1;
        //                   es-=1;
        //                }
        //             else
        //                while (es > 0)
        //                {
        //                   if (nblock >= nblockMAX)
        //                      RETURN(BZ_DATA_ERROR);
        //                   s.tt[nblock] = (UInt32)uc;
        //                   nblock+=1;
        //                   es-=1;
        //                };

        //             continue;
        //          }
        //          else
        //          {

        //             if (nblock >= nblockMAX)
        //                RETURN(BZ_DATA_ERROR);

        //             /*-=1 uc = MTF ( nextSym-1 ) -=1*/
        //             {
        //                Int32 ii, jj, kk, pp, lno, off;
        //                UInt32 nn;
        //                nn = (UInt32)(nextSym - 1);

        //                if (nn < MTFL_SIZE)
        //                {
        //                   /* avoid general-case expense */
        //                   pp = s.mtfbase[0];
        //                   uc = s.mtfa[pp + nn];
        //                   while (nn > 3)
        //                   {
        //                      Int32 z = pp + nn;
        //                      s.mtfa[(z)] = s.mtfa[(z)-1];
        //                      s.mtfa[(z)-1] = s.mtfa[(z)-2];
        //                      s.mtfa[(z)-2] = s.mtfa[(z)-3];
        //                      s.mtfa[(z)-3] = s.mtfa[(z)-4];
        //                      nn -= 4;
        //                   }
        //                   while (nn > 0)
        //                   {
        //                      s.mtfa[(pp + nn)] = s.mtfa[(pp + nn) - 1];
        //                      nn-=1;
        //                   };
        //                   s.mtfa[pp] = uc;
        //                }
        //                else
        //                {
        //                   /* general case */
        //                   lno = nn / MTFL_SIZE;
        //                   off = nn % MTFL_SIZE;
        //                   pp = s.mtfbase[lno] + off;
        //                   uc = s.mtfa[pp];
        //                   while (pp > s.mtfbase[lno])
        //                   {
        //                      s.mtfa[pp] = s.mtfa[pp - 1];
        //                      pp-=1;
        //                   };
        //                   s.mtfbase[lno]+=1;
        //                   while (lno > 0)
        //                   {
        //                      s.mtfbase[lno]-=1;
        //                      s.mtfa[s.mtfbase[lno]] = s.mtfa[s.mtfbase[lno - 1] + MTFL_SIZE - 1];
        //                      lno-=1;
        //                   }
        //                   s.mtfbase[0]-=1;
        //                   s.mtfa[s.mtfbase[0]] = uc;
        //                   if (s.mtfbase[0] == 0)
        //                   {
        //                      kk = MTFA_SIZE - 1;
        //                      for (ii = 256 / MTFL_SIZE - 1; ii >= 0; ii-=1)
        //                      {
        //                         for (jj = MTFL_SIZE - 1; jj >= 0; jj-=1)
        //                         {
        //                            s.mtfa[kk] = s.mtfa[s.mtfbase[ii] + jj];
        //                            kk-=1;
        //                         }
        //                         s.mtfbase[ii] = kk + 1;
        //                      }
        //                   }
        //                }
        //             }
        //             /*-=1 end uc = MTF ( nextSym-1 ) -=1*/
        //             s.unzftab[s.seqToUnseq[uc]]+=1;
        //             if (s.smallDecompress)
        //                s.ll16[nblock] = (UInt16)(s.seqToUnseq[uc]);
        //             else
        //                s.tt[nblock] = (UInt32)(s.seqToUnseq[uc]);
        //             nblock+=1;

        //             GET_MTF_VAL(BZ_X_MTF_5, BZ_X_MTF_6, nextSym);
        //             continue;
        //          }
        //       }

        //       /* Now we know what nblock is, we can do a better sanity
        //          check on s.origPtr.
        //       */
        //       if (s.origPtr < 0 || s.origPtr >= nblock)
        //          RETURN(BZ_DATA_ERROR);

        //       /*-=1 Set up cftab to facilitate generation of T^(-1) -=1*/
        //       /* Check: unzftab entries in range. */
        //       for (i = 0; i <= 255; i+=1)
        //       {
        //          if (s.unzftab[i] < 0 || s.unzftab[i] > nblock)
        //             RETURN(BZ_DATA_ERROR);
        //       }
        //       /* Actually generate cftab. */
        //       s.cftab[0] = 0;
        //       for (i = 1; i <= 256; i+=1)
        //          s.cftab[i] = s.unzftab[i - 1];
        //       for (i = 1; i <= 256; i+=1)
        //          s.cftab[i] += s.cftab[i - 1];
        //       /* Check: cftab entries in range. */
        //       for (i = 0; i <= 256; i+=1)
        //       {
        //          if (s.cftab[i] < 0 || s.cftab[i] > nblock)
        //          {
        //             /* s.cftab[i] can legitimately be == nblock */
        //             RETURN(BZ_DATA_ERROR);
        //          }
        //       }
        //       /* Check: cftab entries non-descending. */
        //       for (i = 1; i <= 256; i+=1)
        //       {
        //          if (s.cftab[i - 1] > s.cftab[i])
        //          {
        //             RETURN(BZ_DATA_ERROR);
        //          }
        //       }

        //       s.state_out_len = 0;
        //       s.state_out_ch = 0;
        //       BZ_INITIALISE_CRC(s.calculatedBlockCRC);
        //       s.state = BZ_X_OUTPUT;
        //       if (s.verbosity >= 2)
        //          VPrintf0("rt+rld");

        //       if (s.smallDecompress)
        //       {

        //          /*-=1 Make a copy of cftab, used in generation of T -=1*/
        //          for (i = 0; i <= 256; i+=1)
        //             s.cftabCopy[i] = s.cftab[i];

        //          /*-=1 compute the T vector -=1*/
        //          for (i = 0; i < nblock; i+=1)
        //          {
        //             uc = (UChar)(s.ll16[i]);
        //             SET_LL(i, s.cftabCopy[uc]);
        //             s.cftabCopy[uc]+=1;
        //          }

        //          /*-=1 Compute T^(-1) by pointer reversal on T -=1*/
        //          i = s.origPtr;
        //          j = GET_LL(i);
        //          do
        //          {
        //             Int32 tmp = GET_LL(j);
        //             SET_LL(j, i);
        //             i = j;
        //             j = tmp;
        //          } while (i != s.origPtr);

        //          s.tPos = s.origPtr;
        //          s.nblock_used = 0;
        //          if (s.blockRandomised)
        //          {
        //             BZ_RAND_INIT_MASK;
        //             BZ_GET_SMALL(s.k0);
        //             s.nblock_used+=1;
        //             BZ_RAND_UPD_MASK;
        //             s.k0 ^= BZ_RAND_MASK;
        //          }
        //          else
        //          {
        //             BZ_GET_SMALL(s.k0);
        //             s.nblock_used+=1;
        //          }
        //       }
        //       else
        //       {

        //          /*-=1 compute the T^(-1) vector -=1*/
        //          for (i = 0; i < nblock; i+=1)
        //          {
        //             uc = (UChar)(s.tt[i] & 0xff);
        //             s.tt[s.cftab[uc]] |= (i << 8);
        //             s.cftab[uc]+=1;
        //          }

        //          s.tPos = s.tt[s.origPtr] >> 8;
        //          s.nblock_used = 0;
        //          if (s.blockRandomised)
        //          {
        //             BZ_RAND_INIT_MASK;
        //             BZ_GET_FAST(s.k0);
        //             s.nblock_used+=1;
        //             BZ_RAND_UPD_MASK;
        //             s.k0 ^= BZ_RAND_MASK;
        //          }
        //          else
        //          {
        //             BZ_GET_FAST(s.k0);
        //             s.nblock_used+=1;
        //          }
        //       }

        //       RETURN(BZ_OK);

        //    endhdr_2:

        //       GET_UCHAR(BZ_X_ENDHDR_2, uc);
        //       if (uc != 0x72)
        //          RETURN(BZ_DATA_ERROR);
        //       GET_UCHAR(BZ_X_ENDHDR_3, uc);
        //       if (uc != 0x45)
        //          RETURN(BZ_DATA_ERROR);
        //       GET_UCHAR(BZ_X_ENDHDR_4, uc);
        //       if (uc != 0x38)
        //          RETURN(BZ_DATA_ERROR);
        //       GET_UCHAR(BZ_X_ENDHDR_5, uc);
        //       if (uc != 0x50)
        //          RETURN(BZ_DATA_ERROR);
        //       GET_UCHAR(BZ_X_ENDHDR_6, uc);
        //       if (uc != 0x90)
        //          RETURN(BZ_DATA_ERROR);

        //       s.storedCombinedCRC = 0;
        //       GET_UCHAR(BZ_X_CCRC_1, uc);
        //       s.storedCombinedCRC = (s.storedCombinedCRC << 8) | ((UInt32)uc);
        //       GET_UCHAR(BZ_X_CCRC_2, uc);
        //       s.storedCombinedCRC = (s.storedCombinedCRC << 8) | ((UInt32)uc);
        //       GET_UCHAR(BZ_X_CCRC_3, uc);
        //       s.storedCombinedCRC = (s.storedCombinedCRC << 8) | ((UInt32)uc);
        //       GET_UCHAR(BZ_X_CCRC_4, uc);
        //       s.storedCombinedCRC = (s.storedCombinedCRC << 8) | ((UInt32)uc);

        //       s.state = BZ_X_IDLE;
        //       RETURN(BZ_STREAM_END);

        //    default:
        //       AssertH(False, 4001);

        break 'outer;
    }
    s.save_i = i;
    s.save_j = j;
    s.save_t = t;
    s.save_alphaSize = alphaSize;
    s.save_nGroups = nGroups;
    s.save_nSelectors = nSelectors;
    s.save_EOB = EOB;
    s.save_groupNo = groupNo;
    s.save_groupPos = groupPos;
    s.save_nextSym = nextSym;
    s.save_nblockMAX = nblockMAX;
    s.save_nblock = nblock;
    s.save_es = es;
    s.save_N = N;
    s.save_curr = curr;
    s.save_zt = zt;
    s.save_zn = zn;
    s.save_zvec = zvec;
    s.save_zj = zj;
    s.save_gSel = gSel;
    s.save_gMinlen = gMinlen;
    s.save_gLimit = gLimit;
    s.save_gBase = gBase;
    s.save_gPerm = gPerm;

    return retVal;
}

#[no_mangle]
pub extern "C" fn dc_undo_mtf_selector_values(s: &mut DState, nGroups: usize, nSelectors: usize) {
    {
        let mut pos = [0_u8; BZ_N_GROUPS as usize];
        for v in 0..nGroups {
            pos[v] = v as u8;
        }

        for i in 0..nSelectors {
            let mut v = s.selectorMtf[i] as usize;
            let tmp = pos[v];
            while v > 0 {
                pos[v] = pos[v - 1];
                v -= 1;
            }
            pos[0] = tmp;
            s.selector[i] = tmp;
        }
    }
}

#[no_mangle]
pub extern "C" fn dc_create_decoding_tables(s: &mut DState, nGroups: usize, alphaSize: usize) {
    for t in 0..nGroups {
        let mut minLen = 32;
        let mut maxLen = 0;
        for i in 0..alphaSize {
            if s.len[t][i] > maxLen {
                maxLen = s.len[t][i];
            }
            if s.len[t][i] < minLen {
                minLen = s.len[t][i];
            }
        }
        bz2_hb_create_decode_tables(
            &mut s.limit[t],
            &mut s.base[t],
            &mut s.perm[t],
            &mut s.len[t],
            minLen as i32,
            maxLen as i32,
            alphaSize as i32,
        );
        s.minLens[t] = minLen as i32;
    }
}

#[no_mangle]
pub extern "C" fn dc_mtf_init(s: &mut DState) {
    let mut kk = MTFA_SIZE as usize - 1;
    for ii in (0..256 / MTFL_SIZE as usize).rev() {
        for jj in (0..MTFL_SIZE as usize).rev() {
            s.mtfa[kk] = (ii * MTFL_SIZE as usize + jj) as u8;
            kk -= 1;
        }
        s.mtfbase[ii] = kk as i32 + 1;
    }
}

#[no_mangle]
pub extern "C" fn mtf_next_sym(s: &mut DState, nextSym: i32) -> u8 {
    let nn = nextSym - 1;

    if nn < MTFL_SIZE as i32 {
        // avoid general-case expense
        mtf_next_sym_opti(s, nn as u32)
    } else {
        // general case
        mtf_next_sym_general(s, nn as u32)
    }
}

fn mtf_next_sym_opti(s: &mut DState, last: u32) -> u8 {
    let mut nn = last as i32;
    let pp = s.mtfbase[0];
    let uc = s.mtfa[(pp + nn) as usize];
    while nn > 3 {
        let z = pp + nn;
        s.mtfa[z as usize] = s.mtfa[(z - 1) as usize];
        s.mtfa[(z - 1) as usize] = s.mtfa[(z - 2) as usize];
        s.mtfa[(z - 2) as usize] = s.mtfa[(z - 3) as usize];
        s.mtfa[(z - 3) as usize] = s.mtfa[(z - 4) as usize];
        nn -= 4;
    }
    while nn > 0 {
        s.mtfa[(pp + nn) as usize] = s.mtfa[(pp + nn - 1) as usize];
        nn -= 1;
    }
    s.mtfa[pp as usize] = uc;
    uc
}

fn mtf_next_sym_general(s: &mut DState, last: u32) -> u8 {
    let nn = last as i32;
    let mut lno = nn / MTFL_SIZE as i32;
    let off = nn % MTFL_SIZE as i32;
    let mut pp = s.mtfbase[lno as usize] + off;
    let uc = s.mtfa[pp as usize];
    while pp > s.mtfbase[lno as usize] {
        s.mtfa[pp as usize] = s.mtfa[(pp - 1) as usize];
        pp -= 1;
    }
    s.mtfbase[lno as usize] += 1;
    while lno > 0 {
        s.mtfbase[lno as usize] -= 1;
        s.mtfa[s.mtfbase[lno as usize] as usize] =
            s.mtfa[(s.mtfbase[(lno - 1) as usize] + MTFL_SIZE as i32 - 1) as usize];
        lno -= 1;
    }
    s.mtfbase[0] -= 1;
    s.mtfa[s.mtfbase[0] as usize] = uc;
    if s.mtfbase[0] == 0 {
        // let mut kk = (MTFA_SIZE - 1) as i32;
        // for ii in (0..256 / MTFL_SIZE).rev() {
        //     for jj in (0..MTFL_SIZE).rev() {
        //         s.mtfa[kk as usize] = s.mtfa[(s.mtfbase[ii as usize] + jj as i32) as usize];
        //         kk -= 1;
        //     }
        //     s.mtfbase[ii as usize] = kk + 1;
        // }
        let mut kk = (MTFA_SIZE - 1) as i32;
        let mut ii = (256 / MTFL_SIZE - 1) as i32;
        while ii >= 0 {
            let mut jj = (MTFL_SIZE - 1) as i32;
            while jj >= 0 {
                s.mtfa[kk as usize] = s.mtfa[(s.mtfbase[ii as usize] + jj) as usize];
                kk -= 1;
                jj -= 1;
            }
            s.mtfbase[ii as usize] = kk + 1;
            ii -= 1;
        }
    }
    uc
}

#[no_mangle]
pub extern "C" fn dc_smallDecompress(s: &mut DState, nblock: i32) {
    let ll4 = unsafe { from_raw_parts_mut(s.ll4, nblock as usize) };
    let ll16 = unsafe { from_raw_parts_mut(s.ll16, nblock as usize) };

    // Make a copy of cftab, used in generation of T
    for i in 0..257 {
        s.cftabCopy[i] = s.cftab[i];
    }

    // compute the T vector
    for i in 0..nblock as usize {
        let uc = ll16[i] as u8;
        SET_LL(ll16, ll4, i as u8, s.cftabCopy[uc as usize] as u32);
        s.cftabCopy[uc as usize] += 1;
    }

    // Compute T^(-1) by pointer reversal on T
    let mut i = s.origPtr as u32;
    let mut j = GET_LL(ll16, ll4, i as u8);
    loop {
        let tmp = GET_LL(ll16, ll4, j as u8);
        SET_LL(ll16, ll4, j as u8, i);
        i = j;
        j = tmp;
        if i == s.origPtr as u32 {
            break;
        }
    }

    s.tPos = s.origPtr as u32;
    s.nblock_used = 0;
    if s.blockRandomised > 0 {
        // BZ_RAND_INIT_MASK;
        s.rNToGo = 0;
        s.rTPos = 0;
        s.k0 = BZ_GET_SMALL(s, nblock);
        s.nblock_used += 1;
        BZ_RAND_UPD_MASK(s);
        s.k0 ^= if s.rNToGo == 1 { 1 } else { 0 };
    } else {
        s.k0 = BZ_GET_SMALL(s, nblock);
        s.nblock_used += 1;
    }
}

fn SET_LL(ll16: &mut [u16], ll4: &mut [u8], i: u8, n: u32) {
    ll16[i as usize] = n as u16 & 0x0000ffff;
    SET_LL4(ll4, i, n >> 16);
}

fn GET_LL(ll16: &mut [u16], ll4: &mut [u8], i: u8) -> u32 {
    (ll16[i as usize] as u32) | (GET_LL4(ll4, i) << 16)
}

fn SET_LL4(ll4: &mut [u8], i: u8, n: u32) {
    if (i & 0x1) == 0 {
        ll4[(i >> 1) as usize] = (ll4[(i >> 1) as usize] & 0xf0) | n as u8;
    } else {
        ll4[(i >> 1) as usize] = (ll4[(i >> 1) as usize] & 0x0f) | (n << 4) as u8;
    }
}

fn GET_LL4(ll4: &mut [u8], i: u8) -> u32 {
    ((ll4[(i >> 1) as usize] as u32) >> (i << 2) & 0x4) & 0xF
}

fn BZ_GET_SMALL(s: &mut DState, nblock: i32) -> i32 {
    let ll4 = unsafe { from_raw_parts_mut(s.ll4, nblock as usize) };
    let ll16 = unsafe { from_raw_parts_mut(s.ll16, nblock as usize) };

    // c_tPos is unsigned, hence test < 0 is pointless.
    if s.tPos >= 100000 * s.blockSize100k as u32 {
        return 1;
    }
    let result = BZ2_indexIntoF(s.tPos as i32, s.cftab.as_mut_ptr());
    s.tPos = GET_LL(ll16, ll4, s.tPos as u8);
    result
}

fn BZ_RAND_UPD_MASK(s: &mut DState) {
    if s.rNToGo == 0 {
        s.rNToGo = BZ2_rNums[s.rTPos as usize];
        s.rTPos += 1;
        if s.rTPos == 512 {
            s.rTPos = 0;
        }
    }
    s.rNToGo -= 1;
}
