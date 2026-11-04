#!/usr/bin/env python3
"""Exact ten-pair exclusion in the nontrivial conductor-83 stratum.

For each base mask in ``{0,2,3,4,5,8,9,10,11,12}``, this certificate
reconstructs the canonical v45 representative and its bit-16 reflection,
applies an integral determinant-one inversion and a frozen saturated-unit
balancing, derives point-specific rational outward bounds, and exhausts the
bounded determinant-one search.  The reflected candidate set is transferred
by the exact involution

    (a,b,c,d) -> (a,-b,-c,d),

whose denominator and numerator squared-norm elements agree at centers
``x`` and ``-x``.  The normalized quadrilateral is nevertheless tested
separately at every one of the twenty centers.

The reflection matrix has determinant -1.  It lies in ``GL_2(O_L)`` but not
in ``GL_2^+(O_L)`` or ``SL_2(O_L)`` and is used only as an enumeration
bijection.  The inversion and every balancing conjugator lie in
``SL_2(O_L)``.  All decisive arithmetic and comparisons are exact.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from math import floor

import sympy as sp

import d7_ambient_points_v45 as ambient
import d7_r83_masks14_30_v47 as pair_core
import residual_orbit_k17_k38_v42 as core


BASE_MASKS = (0, 2, 3, 4, 5, 8, 9, 10, 11, 12)
MASKS = tuple(mask for base in BASE_MASKS for mask in (base, base + 16))
INDEX = 2
LABEL = "R83:nontrivial"
ALL_POSITIVE = (1,) * core.DEGREE
ZERO = core.ZERO


def qs(text: str) -> tuple[Q, ...]:
    """Compact, inspectable exact field-coordinate literal."""
    return tuple(Q(token) for token in text.split())


def qb(text: str) -> tuple[tuple[Q, ...], ...]:
    """Compact, inspectable six-by-six exact bound literal."""
    answer = tuple(qs(row) for row in text.split(";"))
    assert len(answer) == 6 and all(len(row) == core.DEGREE for row in answer)
    return answer


EXPONENTS = {
    0: (0, 0, -1, -1, -2),
    2: (1, 1, 1, 0, 1),
    3: (1, 1, 0, 0, 0),
    4: (0, 0, 0, 0, -1),
    5: (1, -1, 0, -1, -2),
    8: (0, 0, 0, 0, -1),
    9: (1, 0, 0, 0, -1),
    10: (0, 1, 0, 0, -1),
    11: (1, 0, 0, 0, -1),
    12: (0, 0, 0, -1, -2),
}

# Each value is (center, height squared, modulus squared), in the power basis
# 1,T,...,T^5.  The bit-16 partner has the negative center and the same other
# two invariants.
EXPECTED_INVARIANTS = {
    0: (
        qs("-770/83 -8637/166 -176/83 4407/166 -65/166 -334/83"),
        qs("-3191/13778 -40743/27556 7148/6889 10371/6889 -1669/6889 -1990/6889"),
        qs("-642283/6889 -4044696/6889 3240228/6889 4791660/6889 -352136/6889 -762224/6889"),
    ),
    2: (
        qs("6911/166 3738/83 -5205/166 -3645/166 659/166 102/83"),
        qs("113547/27556 147605/27556 -87441/27556 -51767/13778 15011/27556 8495/13778"),
        qs("9816843/6889 9725047/6889 -12425518/6889 -8593105/6889 3510911/6889 1997079/6889"),
    ),
    3: (
        qs("2347/166 3757/83 -659/83 -1911/83 -9/166 167/83"),
        qs("-113/13778 -28757/27556 1753/27556 9709/13778 -129/13778 -3069/27556"),
        qs("-257411/6889 -4494767/6889 -1902134/6889 2099544/6889 988227/6889 -17916/6889"),
    ),
    4: (
        qs("-1503/83 -18695/166 2161/166 14557/166 -343/166 -1336/83"),
        qs("-9758/6889 -259547/27556 10701/6889 88267/13778 -7797/27556 -7106/6889"),
        qs("-4049287/6889 -27103188/6889 3144228/6889 18155959/6889 -217762/6889 -2825382/6889"),
    ),
    5: (
        qs("-1197/83 -7804/83 -8565/166 8602/83 3517/166 -2060/83"),
        qs("-15243/27556 -115845/27556 -62739/27556 63211/13778 26859/27556 -29679/27556"),
        qs("-14564915/6889 -109122203/6889 -60012266/6889 118315665/6889 25635348/6889 -27670247/6889"),
    ),
    8: (
        qs("-603/166 -501/83 2217/166 204/83 -371/83 -102/83"),
        qs("335/27556 -96/6889 1307/13778 -444/6889 -593/27556 139/6889"),
        qs("-114356/6889 -1324162/6889 -2020140/6889 -680217/6889 1102583/6889 534453/6889"),
    ),
    9: (
        qs("-130/83 704/83 2318/83 1596/83 -1057/83 -1679/166"),
        qs("-1265/27556 -13459/27556 -14909/27556 5781/27556 6969/27556 471/27556"),
        qs("-5020873/6889 -42272227/6889 -56295975/6889 4946211/6889 26804832/6889 7188219/6889"),
    ),
    10: (
        qs("1011/166 3107/83 4294/83 -1327/166 -4127/166 -955/166"),
        qs("-473/6889 -8807/13778 -5994/6889 -519/27556 11801/27556 4119/27556"),
        qs("-9286602/6889 -79140774/6889 -114415849/6889 -2171489/6889 54842091/6889 18458848/6889"),
    ),
    11: (
        qs("2069/166 5530/83 2366/83 -10309/166 -1067/83 2125/166"),
        qs("-2059/27556 -21885/27556 -12591/27556 23645/27556 5331/27556 -1389/6889"),
        qs("-2853527/6889 -22017829/6889 -12813971/6889 23059083/6889 5534139/6889 -5240175/6889"),
    ),
    12: (
        qs("-65/83 1929/83 8709/166 881/83 -1317/83 -881/166"),
        qs("6469/27556 49527/27556 78671/27556 19809/27556 -5058/6889 -8203/27556"),
        qs("864249/6889 6980405/6889 11697859/6889 4135058/6889 -1907552/6889 -948191/6889"),
    ),
}

# Rows are denominator, numerator, a, b, c, d bounds, respectively.
EXPECTED_BOUNDS = {
    0: qb("283/100 297/100 51/50 23/10 29/50 117/100;721/100 619/100 9/50 14/5 24/25 223/100;2257/100 253/50 209/100 311/100 1383/100 616/25;97413/100 18569/100 18/5 1468/25 1583/50 6923/50;221/25 243/100 294/25 64/25 837/100 129/10;38159/100 4451/50 507/25 1206/25 479/25 1449/20"),
    2: qb("69/20 3/4 116/25 13/20 43/50 99/50;881/100 157/100 83/100 79/100 141/100 189/50;1849/100 1999/100 23/50 1103/100 467/50 1451/100;23787/20 1173/25 1643/100 827/50 4689/100 2351/10;181/25 959/100 257/100 907/100 113/20 38/5;46589/100 45/2 9279/100 1359/100 1419/50 6151/50"),
    3: qb("291/100 103/100 107/50 207/100 37/50 27/20;743/100 43/20 19/50 63/25 61/50 257/100;548/25 1459/100 99/100 69/20 54/5 2141/100;50161/50 6429/100 757/100 1324/25 4057/100 15939/100;859/100 7 279/50 283/100 653/100 56/5;19649/50 1541/50 171/4 4351/100 491/20 417/5"),
    4: qb("119/50 4/5 9/10 341/100 59/50 39/20;303/50 83/50 4/25 83/20 97/50 371/100;2689/100 473/25 59/25 21/10 679/100 739/50;81759/100 4957/100 159/50 2181/25 1614/25 23083/100;527/50 227/25 1331/100 43/25 411/100 387/50;32027/100 2377/100 448/25 7167/100 3907/100 12079/100"),
    5: qb("113/50 13/10 27/25 43/50 19/25 161/25;144/25 271/100 1/5 26/25 5/4 123/10;707/25 577/50 49/25 167/20 529/50 223/50;7773/10 8129/100 383/100 437/20 4141/100 38253/50;277/25 277/50 276/25 343/50 32/5 117/50;30449/100 3897/100 2159/100 359/20 1253/50 10008/25"),
    8: qb("61/20 28/25 31/25 33/25 143/100 83/50;389/50 233/100 11/50 161/100 47/20 159/50;523/25 673/50 171/100 541/100 28/5 1729/100;1051 6971/100 439/100 3373/100 7817/100 4933/25;41/5 129/20 963/100 89/20 339/100 181/20;4117/10 1671/50 2477/100 2771/100 4731/100 413/4"),
    9: qb("651/100 51/50 6/5 3/4 37/50 153/50;831/50 53/25 11/50 91/100 121/100 146/25;49/5 1479/100 177/100 481/50 1089/100 47/5;11226/5 1586/25 106/25 949/50 4021/100 36293/100;96/25 709/100 249/25 79/10 659/100 123/25;87949/100 1521/50 1197/50 1559/100 1217/50 1899/10"),
    10: qb("789/100 131/100 103/100 103/100 103/100 59/50;503/25 68/25 19/100 5/4 17/10 113/50;809/100 1151/100 207/100 174/25 31/4 608/25;54371/20 1631/20 363/100 656/25 1131/20 14031/100;317/100 138/25 233/20 143/25 469/100 1273/100;106491/100 391/10 1023/50 539/25 1711/50 3671/50"),
    11: qb("63/25 87/100 29/20 24/25 51/50 43/10;643/100 181/100 13/50 29/25 42/25 821/100;2533/100 173/10 73/50 187/25 197/25 669/100;868 5423/100 257/50 122/5 5559/100 25519/50;993/100 83/10 823/100 123/20 477/100 7/2;17001/50 26 1449/50 401/20 841/25 13353/50"),
    12: qb("91/25 417/100 141/100 89/100 7/10 1;927/100 87/10 1/4 27/25 29/25 191/100;439/25 18/5 151/100 201/25 571/50 2883/100;125223/100 26113/100 497/100 227/10 1917/50 2367/20;172/25 173/100 851/100 661/100 691/100 1509/100;49053/100 12519/100 2803/100 373/20 116/5 6193/100"),
}

# a entries, c entries, coupled bottom rows, c!=0 candidates, c=0 candidates.
EXPECTED_COUNTS = {
    0: (24035, 6633, 2023, 8136, 0),
    2: (24137, 6607, 1967, 8334, 0),
    3: (23993, 6597, 1937, 7794, 0),
    4: (24059, 6599, 2015, 8702, 10),
    5: (23959, 6619, 2083, 8878, 0),
    8: (23967, 6613, 2021, 8750, 6),
    9: (24001, 6607, 2085, 9146, 0),
    10: (24039, 6611, 2017, 9436, 6),
    11: (23973, 6617, 2013, 8360, 0),
    12: (24013, 6635, 2015, 8508, 4),
}

EXPECTED_HISTOGRAMS = {
    0: {(0, 0): 1810, (0, 1): 4964, (0, 2): 876, (0, 3): 346, (0, 4): 20, (0, 5): 14, (1, 0): 78, (1, 1): 12, (1, 2): 12, (1, 4): 4},
    16: {(0, 0): 1736, (0, 1): 5054, (0, 2): 886, (0, 3): 346, (0, 4): 24, (0, 5): 10, (1, 0): 50, (1, 1): 18, (1, 2): 12},
    2: {(0, 0): 1658, (0, 1): 5194, (0, 2): 996, (0, 3): 360, (0, 4): 18, (0, 5): 12, (1, 0): 66, (1, 1): 14, (1, 2): 6, (1, 3): 8, (1, 4): 2},
    18: {(0, 0): 1694, (0, 1): 5220, (0, 2): 1006, (0, 3): 296, (0, 4): 28, (0, 5): 10, (1, 0): 50, (1, 1): 22, (1, 2): 4, (1, 3): 2, (1, 4): 2},
    3: {(0, 0): 1644, (0, 1): 4848, (0, 2): 900, (0, 3): 282, (0, 4): 18, (0, 5): 14, (1, 0): 62, (1, 1): 18, (1, 2): 6, (1, 4): 2},
    19: {(0, 0): 1630, (0, 1): 4874, (0, 2): 896, (0, 3): 288, (0, 4): 16, (0, 5): 8, (1, 0): 60, (1, 1): 14, (1, 2): 2, (1, 3): 4, (1, 5): 2},
    4: {(0, 0): 1806, (0, 1): 5434, (0, 2): 928, (0, 3): 394, (0, 4): 38, (0, 5): 10, (1, 0): 68, (1, 1): 26, (1, 2): 2, (1, 3): 2, (1, 4): 4},
    20: {(0, 0): 1808, (0, 1): 5352, (0, 2): 980, (0, 3): 420, (0, 4): 26, (0, 5): 12, (1, 0): 72, (1, 1): 24, (1, 2): 8, (1, 3): 6, (1, 4): 4},
    5: {(0, 0): 1774, (0, 1): 5618, (0, 2): 972, (0, 3): 392, (0, 4): 26, (0, 5): 20, (1, 0): 52, (1, 1): 16, (1, 2): 6, (1, 4): 2},
    21: {(0, 0): 1870, (0, 1): 5484, (0, 2): 1022, (0, 3): 392, (0, 4): 24, (0, 5): 16, (1, 0): 46, (1, 1): 6, (1, 2): 8, (1, 3): 4, (1, 4): 6},
    8: {(0, 0): 1838, (0, 1): 5528, (0, 2): 886, (0, 3): 386, (0, 4): 28, (0, 5): 8, (1, 0): 54, (1, 1): 18, (1, 2): 2, (1, 4): 6, (1, 5): 2},
    24: {(0, 0): 1798, (0, 1): 5550, (0, 2): 900, (0, 3): 402, (0, 4): 24, (0, 5): 10, (1, 0): 50, (1, 1): 20, (1, 4): 2},
    9: {(0, 0): 1752, (0, 1): 5858, (0, 2): 1042, (0, 3): 348, (0, 4): 26, (0, 5): 14, (1, 0): 82, (1, 1): 14, (1, 2): 2, (1, 3): 8},
    25: {(0, 0): 1806, (0, 1): 5822, (0, 2): 1052, (0, 3): 328, (0, 4): 14, (0, 5): 16, (1, 0): 76, (1, 1): 22, (1, 2): 4, (1, 3): 4, (1, 4): 2},
    10: {(0, 0): 1964, (0, 1): 5880, (0, 2): 1066, (0, 3): 404, (0, 4): 36, (0, 5): 4, (1, 0): 66, (1, 1): 16, (1, 2): 4, (1, 3): 2},
    26: {(0, 0): 1970, (0, 1): 5884, (0, 2): 1046, (0, 3): 392, (0, 4): 36, (0, 5): 16, (1, 0): 76, (1, 1): 12, (1, 2): 2, (1, 3): 4, (1, 4): 4},
    11: {(0, 0): 1658, (0, 1): 5182, (0, 2): 986, (0, 3): 408, (0, 4): 24, (0, 5): 12, (1, 0): 66, (1, 1): 18, (1, 2): 2, (1, 3): 2, (1, 4): 2},
    27: {(0, 0): 1792, (0, 1): 5110, (0, 2): 944, (0, 3): 402, (0, 4): 26, (0, 5): 8, (1, 0): 62, (1, 1): 10, (1, 2): 2, (1, 3): 4},
    12: {(0, 0): 1760, (0, 1): 5214, (0, 2): 1012, (0, 3): 364, (0, 4): 34, (0, 5): 20, (1, 0): 82, (1, 1): 18, (1, 2): 2, (1, 4): 6},
    28: {(0, 0): 1852, (0, 1): 5196, (0, 2): 982, (0, 3): 364, (0, 4): 34, (0, 5): 10, (1, 0): 56, (1, 1): 12, (1, 2): 2, (1, 4): 4},
}

EXPECTED_DIGESTS = {
    0: "45804d2344a890e78e1346763d937d6ad1432cc50ca4490b1e502406685f065d",
    16: "bf73238c70cdf49152d44e9fee300041080e96be78c598d45cdf9b38c55cfe92",
    2: "b24eee6e97452507710752855ba850802909c6bb07fce8590ef91e39f354fe32",
    18: "af37bb608248a19682e7a348b2c8cc20dd861f014138877f21daff0b5d814628",
    3: "de62d03cd5d4c5f403d2495661f8229545070a104bb8ceeb46acd7f77f1b8fc2",
    19: "b80ae682698dd57828230fee4c6b82cc75615fc5e8eda375864ef2840b7b519a",
    4: "46710d59df6c922304a528bba2f38e9affc857e8ad436294765defff0635a4e0",
    20: "2101d71d1d4f392c6bf360651f092e2e90b2a201633506b22816b05cf049c822",
    5: "491b99f68f4b8cb95ea9650a355e737dfb67934f505388c62a7b36e3cf44938e",
    21: "484642f37fbccdbd1ca08e9f090d455e0a4c7417e1fdd8b0d4c7ca57bf0a966e",
    8: "dcc56539ddcc7ca63d1b52401bac336528397525b66981519a887cab063c5d22",
    24: "8ce35e562295148abab3aae0025ff13f118608542097f877adcf7961c5bcc042",
    9: "fddabcc85df1115417de6efe4908b691c695eb06809a1af4d6526d24b7bc05f5",
    25: "99fcae815bea2acc753ed956237497a9e8f3a83e308bfacd732b83ef78a7146d",
    10: "5465692d11f69ff0bf3848eba868528787ef10a74c48aac38fd0b616bd0503d6",
    26: "a00e7cc9b32857d6731c1cd177926ed6f46ad109212a67d7c7cdf99642d87b04",
    11: "81454a94c45e6a015f79431b32d2947cc08f2f1d68107a3af47831f2f131c55d",
    27: "06feaeb73141073f4483109e0f3e385fb96853cc27fce9cd650a2c3158d6a0c1",
    12: "e87d1de9b7f8a963be479adedb8962472fc969b68bc4c051c4f36457f9d57a0b",
    28: "70c63bd13e29a08d4b5d56ee159955cfc8fe118ea64d76d5ccc79cafcf9f8ea0",
}

EXPECTED_PAIR_AGGREGATES = {
    0: "4ccbb78db598145e486f6ce920ce2b242414c20fd1db41d0a48b24a46756f7e6",
    2: "b1260e07bd41712965921d1e2270784e0e68e3a032ca2eca6ddb475911d393e1",
    3: "3205f62662e42a3c4b267f242533242535e463789ce27f27ae11ee12ac868cb4",
    4: "eb8317988da0d674f56fa88997d3adde2eeaa3749ea2839fa94e6b2cb187d656",
    5: "393f2e41d0fdc003649301b07ad0673c86e7dcf31223facf82d7ad1f93a5781b",
    8: "850196c198f67497ebec722b6296026ae43f926fd793d0a9694ad1c693ba137d",
    9: "cdd94a0e3e1871145f93793cbc0bbf1f6facf3aa96bcd3c02b0fe3df1279a7b3",
    10: "d45e2712f974189abdaf4ca4179c969072e9c16a52fe27ccbea0a02a15ebac34",
    11: "b4af0ce6777f6fe9a3cc2ebba3226fa1a8ec8f8f2ba23322ccf3014ccf8e451d",
    12: "3d655c9b51c04c70655bf9d52b1b23be5780bebe5607e2e5250094c7ec25dd47",
}
EXPECTED_BATCH_DIGEST = "be161e65208adb27deadd69d0d63cbe5a41d500ebd331c02ee7bf22157a47278"
EXPECTED_PAIR_BATCH_DIGEST = "17fcca5c094efb163b7eb39636429b5b1b75a56ca6cc4930816742c59607c046"


def to_core(value) -> tuple[Q, ...]:
    polynomial = sp.Poly(ambient.reduce_field(value), ambient.T, domain=sp.QQ)
    return tuple(
        Q(int(polynomial.nth(index).p), int(polynomial.nth(index).q))
        for index in range(core.DEGREE)
    )


def field_power(value, exponent: int):
    """Compute an exact power in L, including negative unit exponents."""
    if exponent < 0:
        inverse = ambient.reduce_field(1 / value)
        ambient.assert_integral(inverse)
        return field_power(inverse, -exponent)
    answer = sp.Integer(1)
    base = value
    while exponent:
        if exponent & 1:
            answer = ambient.reduce_l(answer * base)
        base = ambient.reduce_l(base * base)
        exponent //= 2
    return answer


def reconstruct_pair(base_mask: int):
    """Reconstruct, invert, and balance one bit-16 reflected pair."""
    reflected_mask = base_mask + 16
    base = {
        label: matrix for label, matrix, _ in ambient.module_templates(INDEX)
    }[LABEL]
    assert ambient.reduce_l(
        ambient.SQUARE_CLASS_UNITS[reflected_mask]
        + ambient.SQUARE_CLASS_UNITS[base_mask]
    ) == 0
    oriented = {
        mask: ambient.conjugate_by_unit(base, ambient.SQUARE_CLASS_UNITS[mask])
        for mask in (base_mask, reflected_mask)
    }
    for matrix in oriented.values():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)

    reflection = sp.diag(1, -1)
    assert reflection.det() == -1
    assert ambient.sign_vector(sp.Integer(-1)) == (-1,) * core.DEGREE
    ambient.assert_integral_matrix(reflection)
    assert reflection * reflection == sp.eye(2)
    assert ambient.matrix_reduce(
        oriented[reflected_mask]
        - reflection * oriented[base_mask] * reflection
    ) == sp.zeros(2)

    raw_invariants = {
        mask: pair_core.v46.fixed_invariants(matrix)
        for mask, matrix in oriented.items()
    }
    raw_modulus_norm = core.norm(to_core(raw_invariants[base_mask][2]))
    assert 0 < raw_modulus_norm < 1

    inversion = sp.Matrix([[0, -1], [1, 0]])
    inversion_inverse = sp.Matrix([[0, 1], [-1, 0]])
    assert inversion.det() == inversion_inverse.det() == 1
    assert inversion * inversion_inverse == sp.eye(2)
    ambient.assert_integral_matrix(inversion)
    ambient.assert_integral_matrix(inversion_inverse)
    inverted = {
        mask: ambient.matrix_multiply(
            ambient.matrix_multiply(inversion, matrix), inversion_inverse
        )
        for mask, matrix in oriented.items()
    }
    for mask, matrix in inverted.items():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)
        raw_center, raw_y2, raw_modulus2 = raw_invariants[mask]
        center, y2, modulus2 = pair_core.v46.fixed_invariants(matrix)
        assert ambient.reduce_field(center + raw_center / raw_modulus2) == 0
        assert ambient.reduce_field(y2 - raw_y2 / raw_modulus2**2) == 0
        assert ambient.reduce_field(modulus2 - 1 / raw_modulus2) == 0
    assert ambient.matrix_reduce(
        inverted[reflected_mask]
        - reflection * inverted[base_mask] * reflection
    ) == sp.zeros(2)

    generators = ambient.XI_UNITS[:4] + (ambient.POSITIVE_NONSQUARE,)
    unit = sp.Integer(1)
    for generator, exponent in zip(generators, EXPONENTS[base_mask]):
        unit = ambient.reduce_l(unit * field_power(generator, exponent))
    unit_inverse = ambient.reduce_field(1 / unit)
    assert abs(ambient.absolute_norm(unit)) == 1
    assert ambient.reduce_l(unit * unit_inverse - 1) == 0
    ambient.assert_integral(unit)
    ambient.assert_integral(unit_inverse)
    conjugator = sp.diag(unit, unit_inverse)
    conjugator_inverse = sp.diag(unit_inverse, unit)
    assert ambient.reduce_l(conjugator.det() - 1) == 0
    assert ambient.reduce_l(conjugator_inverse.det() - 1) == 0
    ambient.assert_integral_matrix(conjugator)
    ambient.assert_integral_matrix(conjugator_inverse)

    balanced = {
        mask: ambient.matrix_multiply(
            ambient.matrix_multiply(conjugator, matrix), conjugator_inverse
        )
        for mask, matrix in inverted.items()
    }
    square = ambient.reduce_l(unit * unit)
    fourth = ambient.reduce_l(square * square)
    invariants = {}
    for mask, matrix in balanced.items():
        ambient.assert_integral_matrix(matrix)
        ambient.assert_characteristic(matrix, INDEX)
        center, y2, modulus2 = pair_core.v46.fixed_invariants(matrix)
        inverted_center, inverted_y2, inverted_modulus2 = (
            pair_core.v46.fixed_invariants(inverted[mask])
        )
        assert ambient.reduce_l(center - square * inverted_center) == 0
        assert ambient.reduce_l(y2 - fourth * inverted_y2) == 0
        assert ambient.reduce_l(modulus2 - fourth * inverted_modulus2) == 0
        invariants[mask] = tuple(to_core(value) for value in (center, y2, modulus2))
        assert abs(ambient.absolute_norm(matrix[1, 0])) == 6889
    assert ambient.matrix_reduce(
        balanced[reflected_mask]
        - reflection * balanced[base_mask] * reflection
    ) == sp.zeros(2)

    expected_center, expected_y2, expected_modulus2 = EXPECTED_INVARIANTS[base_mask]
    assert invariants[base_mask] == (
        expected_center,
        expected_y2,
        expected_modulus2,
    )
    assert invariants[reflected_mask] == (
        core.neg(expected_center),
        expected_y2,
        expected_modulus2,
    )
    assert core.sign_vector(expected_y2) == ALL_POSITIVE
    assert core.sign_vector(expected_modulus2) == ALL_POSITIVE
    assert core.norm(expected_y2) == Q(1, 28217344)
    assert core.norm(expected_modulus2) == Q(2097130411225081, 6889)
    return balanced, invariants


def outward_hundredth(value: Q) -> Q:
    """Return the next hundredth strictly outside an exact upper endpoint."""
    answer = Q(floor(100 * value) + 1, 100)
    assert answer > value
    return answer


def derive_bounds(base_mask: int, y2, modulus2):
    """Derive and freeze all point-specific bounds by rational intervals."""
    bounds = [[] for _ in range(6)]
    for index in range(core.DEGREE):
        y_lo, y_hi = core.sqrt_interval(
            core.interval_eval(y2, core.ROOT_INTERVALS[index])
        )
        modulus_lo, modulus_hi = core.sqrt_interval(
            core.interval_eval(modulus2, core.ROOT_INTERVALS[index])
        )
        assert y_lo > 0 and modulus_lo > 0
        denominator = core.sqrt_interval(
            (y_hi / core.M_LOWER[index], y_hi / core.M_LOWER[index])
        )[1]
        c_entry = core.sqrt_interval(
            (
                1 / (core.M_LOWER[index] * y_lo),
                1 / (core.M_LOWER[index] * y_lo),
            )
        )[1]
        numerator = core.U_UPPER[index] * denominator
        a_entry = numerator / y_lo
        d_entry = denominator + c_entry * modulus_hi
        b_entry = numerator + a_entry * modulus_hi
        for target, value in zip(
            bounds, (denominator, numerator, a_entry, b_entry, c_entry, d_entry)
        ):
            target.append(outward_hundredth(value))
    answer = tuple(tuple(row) for row in bounds)
    assert answer == EXPECTED_BOUNDS[base_mask]
    return answer


def within_integral_box_exact(x, bounds, equality_hits: set | None = None):
    """Certified compact-box test with an equality-safe Sturm fallback.

    The fixed-scale filter is decisive except when its outward interval meets
    an endpoint.  In that case ``B+x`` and ``B-x`` are exact elements of L;
    their signs at each isolated real root decide membership.  Equality is
    accepted because the search box is compact.
    """
    try:
        return core.within_integral_box(x, bounds)
    except ArithmeticError as error:
        payload = error.args[0] if error.args else ()
        assert (
            isinstance(payload, tuple)
            and payload
            and payload[0]
            in {"left box boundary unresolved", "right box boundary unresolved"}
        ), error

    value = core.as_field(x)
    for index, bound in enumerate(bounds):
        bound_element = core.scale(bound, core.ONE)
        left = core.add(bound_element, value)
        right = core.sub(bound_element, value)
        left_sign = core.sign_at(left, index)
        right_sign = core.sign_at(right, index)
        if left_sign < 0 or right_sign < 0:
            return False
        if equality_hits is not None and left_sign == 0:
            equality_hits.add(("left", x, index, bound))
        if equality_hits is not None and right_sign == 0:
            equality_hits.add(("right", x, index, bound))
    return True


def entry_lists(base_mask: int, y2, bounds):
    """Exhaust the two exact integral entry boxes."""
    pair_core.install_point(y2, bounds)
    equality_hits: set = set()
    a_values = core.enumerate_trace_ellipsoid(
        ZERO,
        sum((value * value for value in bounds[2]), Q(0)),
        lambda value: within_integral_box_exact(value, bounds[2], equality_hits),
    )
    c_values = core.enumerate_trace_ellipsoid(
        ZERO,
        sum((value * value for value in bounds[4]), Q(0)),
        lambda value: within_integral_box_exact(value, bounds[4], equality_hits),
    )
    expected_a, expected_c, _, _, _ = EXPECTED_COUNTS[base_mask]
    assert (len(a_values), len(c_values)) == (expected_a, expected_c)
    assert {core.neg(value) for value in a_values} == set(a_values)
    assert {core.neg(value) for value in c_values} == set(c_values)
    mask3_hits = {
        ("left", (-7, 0, 0, 0, 0, 0), 1, Q(7)),
        ("right", (7, 0, 0, 0, 0, 0), 1, Q(7)),
    }
    if base_mask == 3:
        assert equality_hits == mask3_hits
    else:
        assert not equality_hits
    return a_values, c_values


def bounded_candidates(base_mask: int, x, a_values, c_values):
    """Enumerate the complete bounded determinant-one set, including c=0."""
    den_squares = core.square_bounds(core.DEN_BOUND)
    num_squares = core.square_bounds(core.NUM_BOUND)
    den_radius = sum(den_squares, Q(0))
    num_radius = sum(num_squares, Q(0))
    bottom_rows = []
    x_boxes = tuple(core.scaled_field_box(x, index) for index in range(core.DEGREE))
    y2_boxes = tuple(
        core.scaled_field_box(core.Y2, index) for index in range(core.DEGREE)
    )
    a_boxes_cache = {
        a: tuple(core.scaled_eval_integral(a, index) for index in range(core.DEGREE))
        for a in a_values
    }

    def numerator_bounded(a, b):
        unresolved = False
        a_boxes = a_boxes_cache.get(a)
        if a_boxes is None:
            a_boxes = tuple(
                core.scaled_eval_integral(a, index)
                for index in range(core.DEGREE)
            )
        for index in range(core.DEGREE):
            b_box = core.scaled_eval_integral(b, index)
            real_box = core.scaled_add(
                core.scaled_mul(a_boxes[index], x_boxes[index]), b_box
            )
            square_box = core.scaled_add(
                core.scaled_mul(real_box, real_box),
                core.scaled_mul(
                    core.scaled_mul(a_boxes[index], a_boxes[index]),
                    y2_boxes[index],
                ),
            )
            decision = core.scaled_le(square_box, num_squares[index])
            if decision is False:
                return False
            unresolved |= decision is None
        if not unresolved:
            return True
        return core.all_le_bounds(
            core.numerator_square(core.as_field(a), core.as_field(b), x),
            num_squares,
        )

    for c_index, c in enumerate(c_values):
        cq = core.as_field(c)
        shifted_center = core.neg(core.mul(cq, x))
        c_boxes = tuple(
            core.scaled_eval_integral(c, index) for index in range(core.DEGREE)
        )

        def denominator_predicate(d):
            unresolved = False
            for index in range(core.DEGREE):
                d_box = core.scaled_eval_integral(d, index)
                real_box = core.scaled_add(
                    core.scaled_mul(c_boxes[index], x_boxes[index]), d_box
                )
                square_box = core.scaled_add(
                    core.scaled_mul(real_box, real_box),
                    core.scaled_mul(
                        core.scaled_mul(c_boxes[index], c_boxes[index]),
                        y2_boxes[index],
                    ),
                )
                decision = core.scaled_le(square_box, den_squares[index])
                if decision is False:
                    return False
                unresolved |= decision is None
            if not unresolved:
                return True
            return core.all_le_bounds(
                core.denominator_square(cq, core.as_field(d), x), den_squares
            )

        bottom_rows.extend(
            (c, d)
            for d in core.enumerate_trace_ellipsoid(
                shifted_center, den_radius, denominator_predicate
            )
        )
        if c_index and c_index % 1000 == 0:
            print(
                f"  pair {base_mask}/{base_mask + 16}: "
                f"c entries {c_index}/{len(c_values)}; "
                f"bottom rows={len(bottom_rows)}",
                flush=True,
            )

    nonzero = []
    nonzero_rows = [(c, d) for c, d in bottom_rows if c != ZERO]
    for row_index, (c, d) in enumerate(nonzero_rows):
        cq, dq = core.as_field(c), core.as_field(d)
        affine_data = core.affine_b_data(c, d)
        for a in a_values:
            b = core.affine_integral_b_fast(a, affine_data)
            if b is None or not numerator_bounded(a, b):
                continue
            aq, bq = core.as_field(a), core.as_field(b)
            assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
            nonzero.append((a, b, c, d))
        if row_index and row_index % 500 == 0:
            print(
                f"  pair {base_mask}/{base_mask + 16}: determinant rows "
                f"{row_index}/{len(nonzero_rows)}; candidates={len(nonzero)}",
                flush=True,
            )

    zero = []
    for c, d in bottom_rows:
        if c != ZERO or d == ZERO:
            continue
        dq = core.as_field(d)
        aq = core.inverse(dq)
        if not core.integral(aq):
            continue
        a = core.integral_tuple(aq)
        if not within_integral_box_exact(a, core.A_BOUND):
            continue
        shifted_center = core.neg(core.mul(aq, x))
        for b in core.enumerate_trace_ellipsoid(
            shifted_center,
            num_radius,
            lambda candidate: numerator_bounded(a, candidate),
        ):
            assert core.mul(aq, dq) == core.ONE
            zero.append((a, b, c, d))

    _, _, expected_bottom, expected_nonzero, expected_zero = EXPECTED_COUNTS[base_mask]
    assert (len(bottom_rows), len(nonzero), len(zero)) == (
        expected_bottom,
        expected_nonzero,
        expected_zero,
    )
    assert all(gamma[2] != ZERO for gamma in nonzero)
    assert all(gamma[2] == ZERO for gamma in zero)
    candidates = nonzero + zero
    assert len(set(candidates)) == len(candidates)
    return candidates


def reflect_raw(gamma):
    a, b, c, d = gamma
    return a, tuple(-value for value in b), tuple(-value for value in c), d


def exact_square(coefficient, constant, x, y2):
    """Return (coefficient*x+constant)^2 + coefficient^2*y2 in L."""
    coefficient_q = core.as_field(coefficient)
    constant_q = core.as_field(constant)
    real = core.add(core.mul(coefficient_q, x), constant_q)
    return core.add(
        core.mul(real, real),
        core.mul(core.mul(coefficient_q, coefficient_q), y2),
    )


def reflected_candidates(base_mask: int, x, y2, candidates):
    """Certify the reflection involution, converse, and norm identities."""
    reflected = [reflect_raw(gamma) for gamma in candidates]
    assert len(set(reflected)) == len(reflected) == len(candidates)
    assert [reflect_raw(gamma) for gamma in reflected] == candidates
    for gamma, partner in zip(candidates, reflected):
        a, b, c, d = gamma
        ar, br, cr, dr = partner
        aq, bq, cq, dq = map(core.as_field, gamma)
        arq, brq, crq, drq = map(core.as_field, partner)
        assert core.sub(core.mul(aq, dq), core.mul(bq, cq)) == core.ONE
        assert core.sub(core.mul(arq, drq), core.mul(brq, crq)) == core.ONE
        # These exact field equalities are stronger than equality of their
        # six real norm coordinates and prove both directions preserve every
        # coupled denominator/numerator bound.
        assert exact_square(c, d, x, y2) == exact_square(
            cr, dr, core.neg(x), y2
        )
        assert exact_square(a, b, x, y2) == exact_square(
            ar, br, core.neg(x), y2
        )
        assert (c == ZERO) == (cr == ZERO)
    expected_total = EXPECTED_COUNTS[base_mask][3] + EXPECTED_COUNTS[base_mask][4]
    assert len(reflected) == expected_total
    return reflected


def classify(mask: int, x, candidates):
    """Test the normalized quadrilateral exactly at one geometric center."""
    exclusions = []
    survivors = []
    histogram = {}
    for gamma in candidates:
        inside, failure = core.in_all_coordinate_polygons(gamma, x)
        if inside:
            survivors.append(gamma)
        else:
            assert failure is not None
            histogram[failure] = histogram.get(failure, 0) + 1
            exclusions.append((gamma, failure))
    assert survivors == []
    assert len(exclusions) == len(candidates)
    assert histogram == EXPECTED_HISTOGRAMS[mask]
    assert sum(histogram.values()) == len(candidates)
    digest = sha256()
    for gamma, failure in sorted(exclusions):
        flat = tuple(coefficient for entry in gamma for coefficient in entry)
        digest.update(
            (",".join(map(str, flat)) + f";{failure[0]},{failure[1]}\n").encode(
                "ascii"
            )
        )
    answer = digest.hexdigest()
    assert answer == EXPECTED_DIGESTS[mask]
    return histogram, answer


def check_fibers() -> None:
    """Freeze xor-48 root fibers and their mutual/prior disjointness."""
    assert ambient.PAIR_XOR[LABEL] == 48
    assert set(MASKS) == set(EXPECTED_HISTOGRAMS) == set(EXPECTED_DIGESTS)
    assert all((LABEL, mask) in ambient.POINT_RECORDS for mask in MASKS)
    new_fibers = tuple(frozenset((mask, mask ^ 48)) for mask in MASKS)
    assert len(set(new_fibers)) == len(new_fibers) == 20
    for index, fiber in enumerate(new_fibers):
        assert all(fiber.isdisjoint(other) for other in new_fibers[index + 1 :])
    prior_fibers = tuple(
        frozenset(fiber)
        for fiber in (
            {14, 62}, {30, 46}, {15, 63}, {31, 47},
            {6, 54}, {22, 38}, {1, 49}, {17, 33},
        )
    )
    assert all(
        new_fiber.isdisjoint(prior_fiber)
        for new_fiber in new_fibers
        for prior_fiber in prior_fibers
    )


def main() -> None:
    print("EXACT NONTRIVIAL-R83 TEN-PAIR ORBIT CERTIFICATE v51")
    check_fibers()
    observed_digests = {}
    observed_pair_aggregates = {}
    for pair_index, base_mask in enumerate(BASE_MASKS, start=1):
        reflected_mask = base_mask + 16
        print(
            f"PAIR {pair_index}/10: masks {base_mask}/{reflected_mask}; "
            f"balancing exponents={EXPONENTS[base_mask]}",
            flush=True,
        )
        _, invariants = reconstruct_pair(base_mask)
        x, y2, modulus2 = invariants[base_mask]
        reflected_x, reflected_y2, reflected_modulus2 = invariants[reflected_mask]
        assert reflected_x == core.neg(x)
        assert reflected_y2 == y2 and reflected_modulus2 == modulus2
        bounds = derive_bounds(base_mask, y2, modulus2)
        a_values, c_values = entry_lists(base_mask, y2, bounds)
        print(
            f"  exact a/c entries={len(a_values)}/{len(c_values)}",
            flush=True,
        )
        candidates = bounded_candidates(base_mask, x, a_values, c_values)
        reflected = reflected_candidates(base_mask, x, y2, candidates)
        base_result = classify(base_mask, x, candidates)
        reflected_result = classify(reflected_mask, reflected_x, reflected)
        observed_digests[base_mask] = base_result[1]
        observed_digests[reflected_mask] = reflected_result[1]
        _, _, bottom, nonzero, zero = EXPECTED_COUNTS[base_mask]
        for mask, result in (
            (base_mask, base_result),
            (reflected_mask, reflected_result),
        ):
            print(
                f"  mask={mask}: bottom={bottom}; c!=0={nonzero}; "
                f"c=0={zero}; survivors=0",
                flush=True,
            )
            print(f"  mask={mask}: histogram={result[0]}", flush=True)
            print(f"  mask={mask}: exclusion digest={result[1]}", flush=True)
        pair_aggregate = sha256(
            "".join(
                f"{mask}:{observed_digests[mask]}\n"
                for mask in (base_mask, reflected_mask)
            ).encode("ascii")
        ).hexdigest()
        assert pair_aggregate == EXPECTED_PAIR_AGGREGATES[base_mask]
        observed_pair_aggregates[base_mask] = pair_aggregate
        print(f"  pair aggregate={pair_aggregate}", flush=True)

    assert observed_digests == EXPECTED_DIGESTS
    batch_digest = sha256(
        "".join(f"{mask}:{observed_digests[mask]}\n" for mask in MASKS).encode(
            "ascii"
        )
    ).hexdigest()
    assert batch_digest == EXPECTED_BATCH_DIGEST
    pair_batch_digest = sha256(
        "".join(
            f"{mask}/{mask + 16}:{observed_pair_aggregates[mask]}\n"
            for mask in BASE_MASKS
        ).encode("ascii")
    ).hexdigest()
    assert pair_batch_digest == EXPECTED_PAIR_BATCH_DIGEST
    print("twenty-mask labelled batch digest:", batch_digest)
    print("ten-pair labelled batch digest:", pair_batch_digest)
    print("PASS: exact exhaustive exclusion for ten reflected pairs / twenty points")


if __name__ == "__main__":
    main()
