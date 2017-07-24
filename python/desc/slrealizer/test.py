import moment

def test_covariance_matrix(self, lensID=None, debug=False, convolve=False):
    if lensID is None:
        print 'No lens system selected for plotting.'
        return
    import random
    filter = 'y'
    while filter == 'y':
        randomIndex = random.randint(0, 200)
        filter = self.observation[randomIndex][1]
        # Now visualize the lens system at the epoch defined by the 
        print(self.observation[randomIndex])
        print(self.catalog.get_lens(lensID))
        print("*******************************")
        desc.slrealizer.covariance_matrix(self.observation[randomIndex],
                                   self.catalog.get_lens(lensID))
    return
