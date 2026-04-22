import unittest
import cogsworth.pop as pop
import cogsworth.obs.xrbs as xrbs
import astropy.units as u

class Test(unittest.TestCase):
    def test_bad_inputs(self):
        """Ensure failure with bad input"""
        it_broke = False
        try:
            xrbs.get_xray_luminosity(10*u.Msun)
        except AttributeError:
            it_broke = True
        self.assertTrue(it_broke)

    def test_xray_luminosity(self):
        """Test x-ray luminosity calculation from bpp and bcm"""
        p = pop.Population(10, final_kstar1=[13, 14],bcm_timestep_conditions=[['dtp=100000.0']],
                                     use_default_BSE_settings=True)
        p.perform_stellar_evolution()
        bcm = p.bcm.drop_duplicates(subset="bin_num", keep='last')
        it_broke = False

        try:
            xrbs.get_xray_luminosity(bcm=bcm)
        except Exception as e:
            print(e)
            it_broke = True
        self.assertFalse(it_broke)