"""
Example unit tests for SLRealizer package
"""
import unittest
import desc.slrealizer

class SLRealizerTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'

    def tearDown(self):
        pass

    def test_run(self):
        foo = desc.slrealizer.SLRealizer(self.message)
        self.assertEquals(foo.run(), self.message)

    def test_failure(self):
        self.assertRaises(TypeError, desc.slrealizer.SLRealizer)
        foo = desc.slrealizer.SLRealizer(self.message)
        self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()
