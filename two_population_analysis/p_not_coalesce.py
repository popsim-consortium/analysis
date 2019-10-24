import stdpopsim
import msprime
from stdpopsim import homo_sapiens
from stdpopsim import drosophila_melanogaster
import numpy as np



model = homo_sapiens.GutenkunstThreePopOutOfAfrica()

#model = drosophila_melanogaster.LiStephanTwoPopulation()

dd = msprime.DemographyDebugger(**model.asdict())
dd.print_history()



#t = 600000 # divergence time for drosophila

t=140e3/25 #divergence time for humans

p = dd.coalescence_rate_trajectory(steps=[t],num_samples=[0,2,0])

print(p[1])

N = -t/(2*np.log(p[1]))

print(N)
