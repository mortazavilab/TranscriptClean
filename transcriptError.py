
# This class represents a possible sequencing error (ie difference from reference genome in a transcript) and keeps track of whether it is corrected and why.
# A Transcript2 object may posess multiple errors

class TranscriptError:
    def __init__(self, transcriptID, position, errorType, size, corrected, reason):

        self.transcriptID = transcriptID
        self.position = position
        self.errorType = errorType
        self.size = size
        self.corrected = corrected
        self.reason = reason


    #def update(self, corrected, reason):
        # This function updates the 'corrected' and 'reason' fields of the object
    #    self.corrected = corrected
    #    self.reason = reason

    def printable(self):
        return "\t".join([self.transcriptID, self.position, self.errorType, str(self.size), self.corrected, self.reason])
