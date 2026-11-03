"""Independent references and guards for the new hydro benchmarks."""
import unittest
import numpy as np
from hydro_suite import analyze, sod_points


class HydroReferences(unittest.TestCase):
    def test_standard_sod_states(self):
        rho,u,p=sod_points(np.array([-1.,0.,1.]),.2)
        np.testing.assert_allclose(rho[[0,2]],[1,.125])
        np.testing.assert_allclose(u[[0,2]],[0,0])
        np.testing.assert_allclose(p[[0,2]],[1,.1])
        self.assertAlmostEqual(rho[1],.4263194281784952,places=12)
        self.assertAlmostEqual(u[1],.9274526200489499,places=12)
        self.assertAlmostEqual(p[1],.30313017805064685,places=12)

    def acoustic_frames(self,moving):
        n=32;a=(np.arange(n)+.5)/n
        x,y=np.meshgrid(a,a);phase=2*np.pi*(x+y)
        frames=[]
        for time in (0.,.25/np.sqrt(2),1/np.sqrt(2)):
            delta=1e-4*np.sinc(1/n)**2*np.sin(phase-(2*np.pi*np.sqrt(2)*time if moving else 0))
            frames.append(dict(time=time,step=len(frames),rho=(1+delta)[None],p=(1/1.4+delta)[None],
                               vx=(delta/np.sqrt(2))[None],vy=(delta/np.sqrt(2))[None],vz=np.zeros((1,n,n))))
        cfg=dict(gamma="1.4",x_lower="0",x_upper="1",y_lower="0",y_upper="1")
        return cfg,[a,a,np.array([0.])],frames

    def test_moving_linear_acoustic_reference(self):
        cfg,axes,frames=self.acoustic_frames(True)
        _,summary=analyze("acoustic",cfg,axes,frames)
        self.assertLess(abs(summary["phase_error"]),1e-11)
        self.assertLess(summary["rho_L1"],1e-15)

    def test_static_wave_cannot_pass_one_period_check(self):
        cfg,axes,frames=self.acoustic_frames(False)
        with self.assertRaisesRegex(ValueError,"phase propagation"):
            analyze("acoustic",cfg,axes,frames)


if __name__=="__main__":unittest.main()
