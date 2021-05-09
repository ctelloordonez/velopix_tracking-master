class ForwardSearch:

    def __init__(self, event):
        self.hits = sort_by_phi(event.hits) # sort hits by phi

    # method that checks previous tracks
    def solve(self):

        '***********Parameters ***************'
        NumberOfPreviousTracks = 200     # Parameter for number of previous tracks to check
        XminY = 0.01 # accepted deviation between the x and y ratio values of the track and a hit h
        YminZ = 0.01 # accepted deviation between the z and y ratio values of the track and a hit h
        XminZ = 0.01 # accepted deviation between the x and z ratio values of the track and a hit h
        angleDifference = 0.005 # accepted deviation between the polar angle in x,y values of the current track and a hit h
        trackGreaterThan = 1 # we only make the current_track a real track if its length is greater than this value
        beginOrEnd = -1 # set to 0 to compare hit polar angle to polar angle of the first hit of the track, or -1 to compare to last hit of the track
        '*************************************'
        
        tracks = [] # list of tracks found
        current_track = [] # list of hits representing the track currently under consideration
        current_track.append(self.hits[0]) # initialize by considering the first hit
        self.hits = self.hits[1:] # Take the first hit out of consideration for forming a track with itself
        skipped = 0 # variable to keep track of the number of modules skipped