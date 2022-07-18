import unittest
import GDModeling 

class Test(unittest.TestCase):
    def test_delta_h(self):
        self.assertEqual(GDModeling.delta_h(0.31350340690241546, 5, 0), 0)
        self.assertEqual(GDModeling.delta_h(0.31350340690241546, 5, -32), 2.191109350089283)
        self.assertEqual(GDModeling.delta_h(0.31350340690241546, 5, -3.4), 0.4915299791716679)
        self.assertNotEqual(GDModeling.delta_h(0.31350340690241546, 5, -3.4), 0.49152)
        self.assertEqual(GDModeling.delta_h(0.31350340690241546, 5, 124.21), 5.411717376441741)
        self.assertNotEqual(GDModeling.delta_h(0.31350340690241546, 5, 124.21), 5.411)
    
    def test_sigma_y(self): # ay = 1.36, by = 0.83
        self.assertEqual(GDModeling.sigma_y(x=0, ay=1.36, by=0.83), 0.0) # 1.36 * pow(abs(0), 0.83)
        self.assertEqual(GDModeling.sigma_y(x=-32, ay=1.36, by=0.83), 24.144231712196305)
        self.assertEqual(GDModeling.sigma_y(x=0, ay=1.36, by=0.83), )
        self.assertEqual(GDModeling.sigma_y(x=0, ay=1.36, by=0.83), )
        self.assertEqual(GDModeling.sigma_y(x=0, ay=1.36, by=0.83), )
        self.assertEqual(GDModeling.sigma_y(x=0, ay=1.36, by=0.83), )

    def test_sigma_y(self): # az = 0.275, by = 0.69
        self.assertEqual(GDModeling.sigma_z(x=0, az=0.275, by=0.69), )
        self.assertEqual(GDModeling.sigma_z(x=-32, az=0.275, by=0.69), )
        self.assertEqual(GDModeling.sigma_z(x=0, az=0.275, by=0.69), )
        self.assertEqual(GDModeling.sigma_z(x=0, az=0.275, by=0.69), )
        self.assertEqual(GDModeling.sigma_z(x=0, az=0.275, by=0.69), )

if __name__ == '__main__':
    unittest.main()

"""
    >>> ((pow(0.31350340690241546, (1/3)))*(pow(abs(-32), (2/3)))*1.6)/5
        2.1911093500892833
    >>> ((pow(0.31350340690241546, (1/3)))*(pow(abs(-3.4), (2/3)))*1.6)/5
        0.4915299791716679
    >>> ((pow(0.31350340690241546, (1/3)))*(pow(abs(124.21), (2/3)))*1.6)/5
        5.411717376441741
"""
