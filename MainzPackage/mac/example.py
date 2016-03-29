import sys
from ROOT import gSystem
gSystem.Load("libMainz_deep_MainzPackage")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing MainzPackage..."

sys.exit(0)

