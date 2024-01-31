import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import TranscriptClean.TranscriptClean as TC
@pytest.mark.unit

class TestCombinedDist(object):

    def test_both_inside(self):
        """ Reference:     ----->|          |<-----
            Transcript:      ----->|      |<-----
            dist_0 = -2, dist_1 = +2, combined dist = 4
        """

        assert TC.combinedJunctionDist(-2, 2) == 4

    def test_left_same_right_inside(self):
        """ Reference:     ----->|          |<-----
            Transcript:    ----->|        |<-----
            dist_0 = 0, dist_1 = +2, combined dist = 2
        """

        assert TC.combinedJunctionDist(0, 2) == 2

    def test_left_outside_right_inside(self):
        """ Reference:     ----->|          |<-----
            Transcript:   ----->|       |<-----
            dist_0 = +1, dist_1 = +4, combined dist = 3
        """

        assert TC.combinedJunctionDist(1, 4) == 3
