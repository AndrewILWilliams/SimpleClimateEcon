# SimpleClimateEcon
Simple cost-benefit optimal model of the coupling between the climate system and an idealized economy

## Example:

```python
import numpy as np
import matplotlib.pyplot as plt
from SimpleClimateEcon import CBA

# Initialise model with default parameters
model = CBA()

# Integrate forward in time
model.forward_integration(nyears=100)

# Plot it
fig, ax = plt.subplots()
model.plot(fig=fig, axs=ax)
```
![Impact of varying the discount rate in the model](Example_figure.pdf)

## Can also do more interesting runs to test sensitivity to parameters such as discount rate
![Impact of varying the discount rate in the model](varying_discount_rate.png)
