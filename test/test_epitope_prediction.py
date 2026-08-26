import unittest
import pandas as pd
from antigenfinder import epitopeprediction


class TestSnpEffParsing(unittest.TestCase):
    def test_epitope_prediction(self):
        df = pd.DataFrame({
            'Peptide': ['MTGGDVAPMGREGVTAMHKL', 'QEAMRRATVEREMELRHKNE', 'EELQVDQLWDALLSRELFRP']
        })

        hla = ['HLA-A*02:01', 'HLA-B*07:02']

        results = epitopeprediction.predict_epitopes(df=df, hla_type=hla)
        print(results)

        self.assertEqual(len(results), 100)

if __name__ == '__main__':
    unittest.main()
