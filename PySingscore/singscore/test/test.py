import unittest, pandas, matplotlib, os
from ..singscore import score,rank, permutate, empiricalpval, \
    plotrankdist, nulldistribution, plotdispersion

from pandas.util.testing import assert_frame_equal

BASE_DIR = os.path.dirname(__file__)or'.'

class SingscoreTestCase(unittest.TestCase):




    def test_score_same_identifiers(self):
        """
        test score() when gene identifiers are the same

        :return: .
        """
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            'test_data/entrez_sample.txt'))

        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            '/test_sigs/tgfb_upDown.txt'))
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header='infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        test_data = os.path.normpath('{}/{}'.format(BASE_DIR,
                                    'output/scored_data.txt'))
        df = pandas.read_csv(test_data, header='infer', sep='\t')
        df = df.rename(columns={'Unnamed: 0': None})
        df = df.set_index(keys=df.columns[0])

        assert_frame_equal((score(up_gene=up, down_gene=down,
                                    sample= sample,
                                  norm_method='theoretical')), df)

    def test_score_upsignatures(self):
        """
        test score() when only up-signatures provided

        :return: .
        """
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                                    'test_data/entrez_sample.txt'))

        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                                   '/test_sigs/tgfb_upDown.txt'))
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header='infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        # get a list of ids
        up = list(up['EntrezID'])
        test_data = os.path.normpath('{}/{}'.format(BASE_DIR,
                                                    'output/scored_data_up.txt'))
        df = pandas.read_csv(test_data, header='infer', sep='\t')
        df = df.rename(columns={'Unnamed: 0': None})
        df = df.set_index(keys=df.columns[0])

        assert_frame_equal((score(up_gene=up,
                                  sample=sample,
                                  norm_method='theoretical')), df)

    def test_rank_same_identifiers(self):
        """
        test rank() when gene identifiers are the same

        :return: .
        """
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            'test_data/entrez_sample.txt'))

        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            '/test_sigs/tgfb_upDown.txt'))
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header='infer', sep='\t')
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        sig.close()
        # data to test df are correct
        test_data = os.path.normpath('{}/{}'.format(BASE_DIR,
                                                    'output/ranked_data.txt'))
        df = pandas.read_csv(test_data, header='infer', sep='\t')
        df = df.rename(columns={'Unnamed: 0': None})
        df = df.set_index(keys=df.columns[0])
        df_columns = df.columns.tolist()

        assert_frame_equal((rank(up_gene=up,
                                 down_gene=down,
                                 sample=sample,
                                 norm_method='theoretical')[df_columns]), df)

    def test_permutation(self):
        """
        test permutations() using only 10 reps for speed.
        :return: .
        """
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            'test_data/entrez_sample.txt'))
        f = open(file_path, 'r')

        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()

        test_data = os.path.normpath('{}/{}'.format(BASE_DIR,
                                                    'output/permd.txt'))
        df = pandas.read_csv(test_data, header='infer', sep='\t')
        df = df.rename(columns={'Unnamed: 0': None})
        df = df.set_index(keys=df.columns[0])
        df.index = range(10)

        assert_frame_equal(permutate(sample=sample, n_up=50, n_down=50,
                                        reps=10),df)


    def test_pvalue(self):
        """
        test empirical_pval() using only 10 reps for speed.
        :return: .
        """
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            'test_data/entrez_sample.txt'))
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            '/test_sigs/tgfb_upDown.txt'))
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header='infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        p = permutate(sample=sample, n_up=50, n_down=50,
                                        reps=10)
        scores = score(up_gene=up, down_gene=down,sample=sample)


        self.assertIsInstance(empiricalpval(permutations=p, score=scores),
                              pandas.DataFrame)

    def test_barcode(self):
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            'test_data/entrez_sample.txt'))
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                                '/test_sigs/tgfb_upDown.txt'))
        sigs = pandas.read_csv(sig_path, header='infer', sep='\t')
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        r = rank(up_gene=up, down_gene=down,
                           sample=sample[['D_Ctrl_R1']],
                           norm_method='theoretical')

        self.assertIsInstance(plotrankdist(ranks=r, show=False),
                              matplotlib.figure.Figure)

    def test_dispersion_plot(self):
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            'test_data/entrez_sample.txt'))
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            '/test_sigs/tgfb_upDown.txt'))
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header='infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])

        scores = score(up_gene=up, down_gene=down,
                                    sample=sample, full_data=True)
        self.assertIsInstance(plotdispersion(score=scores,show=False),
                              matplotlib.figure.Figure)

    def test_null_dist(self):
        file_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            'test_data/entrez_sample.txt'))
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = os.path.normpath('{}/{}'.format(BASE_DIR,
                                            '/test_sigs/tgfb_upDown.txt'))
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header='infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        p = permutate(sample=sample[['D_Ctrl_R1']], n_up=50, n_down=50,
                      reps=10)
        scores = score(up_gene=up, down_gene=down, sample=sample[['D_Ctrl_R1']])

        self.assertIsInstance(nulldistribution(permutations=p, score=scores,
                                               show = False),
                              matplotlib.figure.Figure)


if __name__ == '__main__':
    unittest.main()
