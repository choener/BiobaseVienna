{-# LANGUAGE NoMonomorphismRestriction #-}

module Biobase.Vienna.Export where

import Control.Arrow
import Data.Array.Repa.Index
import Data.Array.Repa.Shape
import Data.List (intersperse)
import Data.List.Split
import qualified Data.Map as M
import Text.Printf
import qualified Data.Vector.Unboxed as VU

import Biobase.Primary
import Biobase.Secondary.Vienna
import Data.PrimitiveArray
import Data.PrimitiveArray.Unboxed.Zero

import Biobase.Vienna



-- * Export as a ViennaRNA 2004 ".par" file

asPar :: Vienna2004 -> Vienna2004 -> String
asPar trnr trnrH = hdr ++ blocks ++ mlps ++ ninios ++ misc ++ triloops ++ tetra ++ hexa ++ "\n#END" where
  hdr = "## RNAfold parameter file v2.0\n\n"
  mlps = printf "# ML_params\n%7d %7d %7d %7d %7d %7d\n\n"
    (multiNuc trnr) (multiNuc trnrH)
    (multiOffset trnr) (multiOffset trnrH)
    (multiHelix trnr) (multiHelix trnrH)
  ninios = printf "# NINIO\n%7d %7d %7d\n\n"
    (ninio trnr) (ninio trnrH) (maxNinio trnr)
  misc = printf "# Misc\n %7d %7d %7d %7d\n\n"
    (intermolecularInit trnr) (intermolecularInit trnrH) (termAU trnr) (termAU trnrH)
  triloops = "# Triloops\n" ++ printHairpinAssocs 5 trnr trnrH ++ "\n"
  tetra = "# Tetraloops\n" ++ printHairpinAssocs 6 trnr trnrH ++ "\n"
  hexa = "# Hexaloops\n" ++ printHairpinAssocs 8 trnr trnrH ++ "\n"
  blocks = concat $ zipWith (++)
    -- entropy terms
    [ printBlock "stack"                7 pp2lkey     $ stack trnr
    , printBlock "mismatch_hairpin"     5 pbb2lkey    $ hairpinMM trnr
    , printBlock "mismatch_interior"    5 pbb2lkey    $ iloopMM trnr
    , printBlock "mismatch_interior_1n" 5 pbb2lkey    $ iloop1xnMM trnr
    , printBlock "mismatch_interior_23" 5 pbb2lkey    $ iloop2x3MM trnr
    , printBlock "mismatch_multi"       5 pbb2lkey    $ multiMM trnr
    , printBlock "mismatch_exterior"    5 pbb2lkey    $ extMM trnr
    , printBlock "dangle5"              5 pb2lkey     $ dangle5 trnr
    , printBlock "dangle3"              5 pb2lkey     $ dangle3 trnr
    , printBlock "int11"                5 ppbb2lkey   $ iloop1x1 trnr
    , printBlock "int21"                5 ppbbb2lkey  $ iloop2x1 trnr
    , printBlock22 "int22"              4 ppbbbb2lkey $ iloop2x2 trnr
    , printLinear "hairpin"            10             $ hairpinL trnr
    , printLinear "bulge"              10             $ bulgeL trnr
    , printLinear "interior"           10             $ iloopL trnr
    ]
    -- enthalpy terms
    [ printBlockH "stack"                7 pp2lkey     $ stack trnrH
    , printBlockH "mismatch_hairpin"     5 pbb2lkey    $ hairpinMM trnrH
    , printBlockH "mismatch_interior"    5 pbb2lkey    $ iloopMM trnrH
    , printBlockH "mismatch_interior_1n" 5 pbb2lkey    $ iloop1xnMM trnrH
    , printBlockH "mismatch_interior_23" 5 pbb2lkey    $ iloop2x3MM trnrH
    , printBlockH "mismatch_multi"       5 pbb2lkey    $ multiMM trnrH
    , printBlockH "mismatch_exterior"    5 pbb2lkey    $ extMM trnrH
    , printBlockH "dangle5"              5 pb2lkey     $ dangle5 trnrH
    , printBlockH "dangle3"              5 pb2lkey     $ dangle3 trnrH
    , printBlockH "int11"                5 ppbb2lkey   $ iloop1x1 trnrH
    , printBlockH "int21"                5 ppbbb2lkey  $ iloop2x1 trnrH
    , printBlock22H "int22"              4 ppbbbb2lkey $ iloop2x2 trnrH
    , printLinearH "hairpin"            10             $ hairpinL trnrH
    , printLinearH "bulge"              10             $ bulgeL trnrH
    , printLinearH "interior"           10             $ iloopL trnrH
    ]



-- * Helper functions

-- | Show the key of the line, minus the changing last key

showKey :: [(LKey,Int)] -> String
showKey xs =
  "   /* " ++
  (concat $ intersperse "," $ init $ (map show ps) ++ (map show ns)) ++
  " */"
  where
    (ps,ns) = fst $ head xs

-- | Transform from tuple-based keys to a pair of list-based keys.

type LKey = ([ViennaPair],[Nuc])

pb2lkey (Z:.p1:.b1) = ([p1],[b1])
pbb2lkey (Z:.p1:.b1:.b2) = ([p1],[b1,b2])
pp2lkey (Z:.k1:.k2) = ([k1,k2],[])
ppbb2lkey (Z:.p1:.p2:.b1:.b2) = ([p1,p2],[b1,b2])
ppbbb2lkey (Z:.p1:.p2:.b1:.b2:.b3) = ([p1,p2],[b1,b2,b3])
ppbbbb2lkey (Z:.p1:.p2:.b1:.b2:.b3:.b4) = ([p1,p2],[b1,b2,b3,b4])

-- | Print a block.

printBlock = printBlockG noNP where
  noNP ((ps,ns),v) = not $ any (==vpNP) ps || any (==nIMI) ns

printBlockH s = printBlock (s ++ "_enthalpies")

printBlock22 = printBlockG noNSNPE where
  noNSNPE ((ps,ns),v) = not $ any (==vpNP) ps || any (==vpNS) ps || any (==nN) ns || any (==nIMI) ns

printBlock22H s = printBlock22 (s ++ "_enthalpies")

printBlockG fltr s k tolkey xs' =
  let
    xs = filter fltr $ map (first tolkey)  $ assocs xs'
  in
    printf "# %s\n" s ++
    (concatMap printLine $ splitEvery k xs) ++
    "\n"

printLine xs =
  concatMap printVal xs ++ " " ++ showKey xs ++
  printf "\n"

printVal (k,v)
  | v > 10000 = printf "   INF"
  | otherwise = printf "%6d" v

-- | A linear block is more boring

printLinear s k xs' = let xs = assocs xs' in
    printf "# %s\n" s ++
    (concatMap (\ys -> concatMap printVal ys ++ "\n") $ splitEvery k xs) ++
    "\n"

printHairpinAssocs l trnr trnrH = res where
  res = concat $ zipWith (\(k,v) vH -> printf "%s %7d %7d\n" (concatMap show k) v vH) xs ys
  xs = filter ((==l).length.fst) $ map (\(k,v) -> (mkString k,v)) $ M.assocs $ hairpinLookup trnr
  ys = map snd $ filter ((==l).length.fst) $ map (\(k,v) -> (mkString k,v)) $ M.assocs $ hairpinLookup trnrH
  mkString = let convT x = if x =='T' then 'U' else x in map (convT . fromNuc) . VU.toList

printLinearH s = printLinear (s ++ "_enthalpies")
