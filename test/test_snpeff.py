import unittest

import antigenfinder.snpeff


class TestSnpEffParsing(unittest.TestCase):
    def test_parsing(self):

        x = antigenfinder.snpeff.parse_consequence('p.Ser81Ser')
        self.assertEqual(1, len(x))
        self.assertEqual(x[0].pos, 81)
        self.assertEqual(x[0].ref, 'S')
        self.assertEqual(x[0].alt, 'S')

        x = antigenfinder.snpeff.parse_consequence('p.Pro491_Gly492del')
        self.assertEqual(2, len(x))
        self.assertEqual(x[0].pos, 491)
        self.assertEqual(x[0].ref, 'P')
        self.assertEqual(x[0].alt, '-')

        self.assertEqual(x[1].pos, 492)
        self.assertEqual(x[1].ref, 'G')
        self.assertEqual(x[1].alt, '-')

        x = antigenfinder.snpeff.parse_consequence('p.Asp1711_Asp1712del')
        self.assertEqual(2, len(x), 'p.Asp1711_Asp1712del')
        self.assertEqual(x[0].pos, 1711, 'p.Asp1711_Asp1712del')
        self.assertEqual(x[0].ref, 'D', 'p.Asp1711_Asp1712del')
        self.assertEqual(x[0].alt, '-', 'p.Asp1711_Asp1712del')

        self.assertEqual(x[1].pos, 1712, 'p.Asp1711_Asp1712del')
        self.assertEqual(x[1].ref, 'D', 'p.Asp1711_Asp1712del')
        self.assertEqual(x[1].alt, '-', 'p.Asp1711_Asp1712del')

        # TODO: see NTs!
        x = antigenfinder.snpeff.parse_consequence('p.Cys48_Gly49insCysSerSerGlyGlyCys')
        self.assertEqual(2, len(x))
        self.assertEqual(x[0].pos, 48)
        self.assertEqual(x[0].ref, 'C')
        self.assertEqual(x[0].alt, 'CSSGGC')

        self.assertEqual(x[1].pos, 49)
        self.assertEqual(x[1].ref, 'G')
        self.assertEqual(x[1].alt, '-')

        x = antigenfinder.snpeff.parse_consequence('p.Trp8_His9delinsCys')
        self.assertEqual(2, len(x))
        self.assertEqual(x[0].pos, 8)
        self.assertEqual(x[0].ref, 'W')
        self.assertEqual(x[0].alt, '-C')

        self.assertEqual(x[1].pos, 9)
        self.assertEqual(x[1].ref, 'H')
        self.assertEqual(x[1].alt, '-')

if __name__ == '__main__':
    unittest.main()
