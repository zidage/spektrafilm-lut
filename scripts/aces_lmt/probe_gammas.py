"""Compare channel printing-density gammas: Status M (profile), APD, and spektrafilm's 2383 print exposure."""
import sys, numpy as np
from spektrafilm_lut_creator.aces_lmt import NegativeADXModel, NegativeADXSpec, ADX_CHANNEL_GAINS
from spektrafilm.runtime.params_builder import init_params, digest_params
from spektrafilm.runtime.pipeline import SimulationPipeline

for film in sys.argv[1:] or ["kodak_vision3_250d"]:
    m = NegativeADXModel(NegativeADXSpec(film_profile=film, balance="dmin"))
    stops = np.array([-2.0, 0.0, 2.0, 4.0, 6.0])
    aces = (0.18 * 2.0 ** stops)[:, None] * np.ones(3)
    cmy = m.cmy_film(aces)
    apd = m.apd_from_cmy(cmy)
    p = init_params(film_profile=film, print_profile="kodak_2383")
    p.debug.lut_mode = True; p.io.input_color_space = "ACES2065-1"
    p = digest_params(p); pipe = SimulationPipeline(p)
    le_print = np.asarray(pipe.process(aces.reshape(1, -1, 3), collect="log_e_print")).reshape(-1, 3)
    st = 0.30103 * 2
    def g(x):  # gamma over [-2,+2] stops and [+2,+6]
        return ((x[2] - x[0]) / (2 * st)).round(3), ((x[4] - x[2]) / (2 * st)).round(3)
    print(film)
    print("  statusM cmy gamma  mid/high", g(cmy))
    print("  APD gamma          mid/high", g(apd), " x ADX gains ->", g(apd * ADX_CHANNEL_GAINS))
    print("  2383 print-density mid/high", g(-le_print))
