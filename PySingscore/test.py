import unittest

if __name__ == '__main__':
    loader = unittest.TestLoader()

    tests = loader.discover(start_dir='./', pattern='test*.py')
    unittest.TextTestRunner().run(tests)