import unittest

from nhe import imt as imt_module


class BaseIMTTestCase(unittest.TestCase):
    class TestIMT(imt_module._IMT):
        _fields = ('foo', 'bar')

        def __new__(cls, foo, bar):
            return imt_module._IMT.__new__(cls, foo, bar)

    def test_base(self):
        self.assertEqual(getattr(self.TestIMT, '__slots__'), ())
        self.assertFalse(hasattr(self.TestIMT(1, 2), '__dict__'))
        imt = self.TestIMT(bar=2, foo=1)
        self.assertEqual(str(imt), 'TestIMT(foo=1, bar=2)')

    def test_equality(self):
        self.assertTrue(self.TestIMT(1, 1) == self.TestIMT(1, 1))
        self.assertTrue(self.TestIMT(1, 5) == self.TestIMT(1, 5))
        self.assertTrue(self.TestIMT('foo', 'bar') ==
                        self.TestIMT('foo', 'bar'))

        self.assertFalse(self.TestIMT(1, 1) != self.TestIMT(1, 1))
        self.assertFalse(self.TestIMT(1, 5) != self.TestIMT(1, 5))
        self.assertFalse(self.TestIMT('foo', 'bar') !=
                         self.TestIMT('foo', 'bar'))

        self.assertFalse(self.TestIMT(1, 1) == (1, 1))
        self.assertFalse(self.TestIMT(0, 1) == self.TestIMT(1, 1))
        self.assertFalse(self.TestIMT(1, 5.4) == self.TestIMT(1, 5))
        self.assertFalse(self.TestIMT('foo', 'bar') ==
                         self.TestIMT('fooz', 'bar'))
        self.assertFalse(self.TestIMT('foo', 'bar') == ('foo', 'bar'))
        self.assertFalse(self.TestIMT(False, False) == False)

        self.assertTrue(self.TestIMT(1, 1) != (1, 1))
        self.assertTrue(self.TestIMT(0, 1) != self.TestIMT(1, 1))
        self.assertTrue(self.TestIMT(1, 5.4) != self.TestIMT(1, 5))
        self.assertTrue(self.TestIMT('foo', 'bar') !=
                        self.TestIMT('fooz', 'bar'))
        self.assertTrue(self.TestIMT('foo', 'bar') != ('foo', 'bar'))
        self.assertTrue(self.TestIMT(False, False) != False)

    def test_hash(self):
        imt1 = self.TestIMT('some', 'thing')
        self.assertEqual(hash(imt1), hash(self.TestIMT('some', 'thing')))
        self.assertNotEqual(hash(imt1), hash(self.TestIMT('other', 'thing')))

        class TestIMT2(self.TestIMT):
            pass

        imt2 = TestIMT2('some', 'thing')
        self.assertNotEqual(hash(imt1), hash(imt2))


class SATestCase(unittest.TestCase):
    def test_wrong_period(self):
        self.assertRaises(ValueError, imt_module.SA, period=0, damping=5)
        self.assertRaises(ValueError, imt_module.SA, period=-1, damping=5)

    def test_wrong_damping(self):
        self.assertRaises(ValueError, imt_module.SA, period=0.1, damping=0)
        self.assertRaises(ValueError, imt_module.SA, period=0.1, damping=-1)
