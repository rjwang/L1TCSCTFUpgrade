# L1TCSCTFUpgrade
==================
General:

	setenv SCRAM_ARCH slc6_amd64_gcc472
        cmsrel CMSSW_6_2_3
        cd CMSSW_6_2_3/src
        cmsenv
        git clone git@github.com:rjwang/L1TCSCTFUpgrade.git ./
        scramv1 b -j 20
