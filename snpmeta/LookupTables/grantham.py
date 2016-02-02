#!/usr/bin/env python
"""Contains a lookup table for amino acid replacement severity from Grantham
1974, Science."""

GSCORES = {
    ('S', 'R'): 110,
    ('S', 'L'): 145,
    ('S', 'P'): 74,
    ('S', 'T'): 58,
    ('S', 'A'): 99,
    ('S', 'V'): 124,
    ('S', 'G'): 56,
    ('S', 'I'): 142,
    ('S', 'F'): 155,
    ('S', 'Y'): 144,
    ('S', 'C'): 112,
    ('S', 'H'): 89,
    ('S', 'Q'): 68,
    ('S', 'N'): 46,
    ('S', 'K'): 121,
    ('S', 'D'): 65,
    ('S', 'E'): 80,
    ('S', 'M'): 135,
    ('S', 'W'): 177,
    ('R', 'L'): 102,
    ('R', 'P'): 103,
    ('R', 'T'): 71,
    ('R', 'A'): 112,
    ('R', 'V'): 96,
    ('R', 'G'): 125,
    ('R', 'I'): 97,
    ('R', 'F'): 97,
    ('R', 'Y'): 77,
    ('R', 'C'): 180,
    ('R', 'H'): 29,
    ('R', 'Q'): 43,
    ('R', 'N'): 86,
    ('R', 'K'): 26,
    ('R', 'D'): 96,
    ('R', 'E'): 54,
    ('R', 'M'): 91,
    ('R', 'W'): 101,
    ('L', 'P'): 98,
    ('L', 'T'): 92,
    ('L', 'A'): 96,
    ('L', 'V'): 32,
    ('L', 'G'): 138,
    ('L', 'I'): 5,
    ('L', 'F'): 22,
    ('L', 'Y'): 36,
    ('L', 'C'): 198,
    ('L', 'H'): 99,
    ('L', 'Q'): 113,
    ('L', 'N'): 153,
    ('L', 'K'): 107,
    ('L', 'D'): 172,
    ('L', 'E'): 138,
    ('L', 'M'): 15,
    ('L', 'W'): 61,
    ('P', 'T'): 38,
    ('P', 'A'): 27,
    ('P', 'V'): 68,
    ('P', 'G'): 42,
    ('P', 'I'): 95,
    ('P', 'F'): 114,
    ('P', 'Y'): 110,
    ('P', 'C'): 169,
    ('P', 'H'): 77,
    ('P', 'Q'): 76,
    ('P', 'N'): 91,
    ('P', 'K'): 103,
    ('P', 'D'): 108,
    ('P', 'E'): 93,
    ('P', 'M'): 87,
    ('P', 'W'): 147,
    ('T', 'A'): 58,
    ('T', 'V'): 69,
    ('T', 'G'): 59,
    ('T', 'I'): 89,
    ('T', 'F'): 103,
    ('T', 'Y'): 92,
    ('T', 'C'): 149,
    ('T', 'H'): 47,
    ('T', 'Q'): 42,
    ('T', 'N'): 65,
    ('T', 'K'): 78,
    ('T', 'D'): 85,
    ('T', 'E'): 65,
    ('T', 'M'): 81,
    ('T', 'W'): 128,
    ('A', 'V'): 64,
    ('A', 'G'): 60,
    ('A', 'I'): 94,
    ('A', 'F'): 113,
    ('A', 'Y'): 112,
    ('A', 'C'): 195,
    ('A', 'H'): 86,
    ('A', 'Q'): 91,
    ('A', 'N'): 111,
    ('A', 'K'): 106,
    ('A', 'D'): 126,
    ('A', 'E'): 107,
    ('A', 'M'): 84,
    ('A', 'W'): 148,
    ('V', 'G'): 109,
    ('V', 'I'): 29,
    ('V', 'F'): 50,
    ('V', 'Y'): 55,
    ('V', 'C'): 192,
    ('V', 'H'): 84,
    ('V', 'Q'): 96,
    ('V', 'N'): 133,
    ('V', 'K'): 97,
    ('V', 'D'): 152,
    ('V', 'E'): 121,
    ('V', 'M'): 21,
    ('V', 'W'): 88,
    ('G', 'I'): 135,
    ('G', 'F'): 153,
    ('G', 'Y'): 147,
    ('G', 'C'): 159,
    ('G', 'H'): 98,
    ('G', 'Q'): 87,
    ('G', 'N'): 80,
    ('G', 'K'): 127,
    ('G', 'D'): 94,
    ('G', 'E'): 98,
    ('G', 'M'): 127,
    ('G', 'W'): 184,
    ('I', 'F'): 21,
    ('I', 'Y'): 33,
    ('I', 'C'): 198,
    ('I', 'H'): 94,
    ('I', 'Q'): 109,
    ('I', 'N'): 149,
    ('I', 'K'): 102,
    ('I', 'D'): 168,
    ('I', 'E'): 134,
    ('I', 'M'): 10,
    ('I', 'W'): 61,
    ('F', 'Y'): 22,
    ('F', 'C'): 205,
    ('F', 'H'): 100,
    ('F', 'Q'): 116,
    ('F', 'N'): 158,
    ('F', 'K'): 102,
    ('F', 'D'): 177,
    ('F', 'E'): 140,
    ('F', 'M'): 28,
    ('F', 'W'): 40,
    ('Y', 'C'): 194,
    ('Y', 'H'): 83,
    ('Y', 'Q'): 99,
    ('Y', 'N'): 143,
    ('Y', 'K'): 85,
    ('Y', 'D'): 160,
    ('Y', 'E'): 122,
    ('Y', 'M'): 36,
    ('Y', 'W'): 37,
    ('C', 'H'): 174,
    ('C', 'Q'): 154,
    ('C', 'N'): 139,
    ('C', 'K'): 202,
    ('C', 'D'): 154,
    ('C', 'E'): 170,
    ('C', 'M'): 196,
    ('C', 'W'): 215,
    ('H', 'Q'): 24,
    ('H', 'N'): 68,
    ('H', 'K'): 32,
    ('H', 'D'): 81,
    ('H', 'E'): 40,
    ('H', 'M'): 87,
    ('H', 'W'): 115,
    ('Q', 'N'): 46,
    ('Q', 'K'): 53,
    ('Q', 'D'): 61,
    ('Q', 'E'): 29,
    ('Q', 'M'): 101,
    ('Q', 'W'): 130,
    ('N', 'K'): 94,
    ('N', 'D'): 23,
    ('N', 'E'): 42,
    ('N', 'M'): 142,
    ('N', 'W'): 174,
    ('K', 'D'): 101,
    ('K', 'E'): 56,
    ('K', 'M'): 95,
    ('K', 'W'): 110,
    ('D', 'E'): 45,
    ('D', 'M'): 160,
    ('D', 'W'): 181,
    ('E', 'M'): 126,
    ('E', 'W'): 152,
    ('M', 'W'): 67,
}
