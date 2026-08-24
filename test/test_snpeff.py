import unittest

import snpeff


class TestSnpEffParsing(unittest.TestCase):
    def test_parsing(self):

        x = snpeff.parse_consequence('p.Ser81Ser')
        self.assertEqual(1, len(x))
        self.assertEqual(x[0].pos_start, 81)
        self.assertEqual(x[0].ref, 'S')
        self.assertEqual(x[0].alt, 'S')


        x = snpeff.parse_consequence('p.Pro491_Gly492del')
        self.assertEqual(1, len(x))
        self.assertEqual(x[0].pos_start, 491)
        self.assertEqual(x[0].ref, 'P')
        self.assertEqual(x[0].alt, '-')

        self.assertEqual(x[1].pos_start, 492)
        self.assertEqual(x[1].ref, 'G')
        self.assertEqual(x[1].alt, '-')

        x = snpeff.parse_consequence('p.Asp1711_Asp1712del')
        print(x)

        x = snpeff.parse_consequence('p.Cys48_Gly49insCysSerSerGlyGlyCys')
        print(x)

        x = snpeff.parse_consequence('p.Trp8_His9delinsCys')
        print(x)


if __name__ == '__main__':
    unittest.main()
