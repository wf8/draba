
	mcmc-iterations	  	500000
	burn-in				50000
	mcmc-sample-skip	10
	iterations-per-log  50
	logs-per-line       10
	start-mig			25000		# 0 is default, need adj after initial runs
	iterations-per-log	1000

	find-finetunes		TRUE
		
	tau-theta-print		10000.0
	tau-theta-alpha		1.0			# for STD/mean ratio of 100%
	tau-theta-beta		10000.0		# for mean of 1e-4

	mig-rate-print		0.001
	mig-rate-alpha		0.002
	mig-rate-beta		0.00001

GENERAL-INFO-END

CURRENT-POPS-START	

	POP-START
		name		Smokie
		samples		DRA_1000 d DRA_1002 d WVA_13 d WVA_12 d
	POP-END
	
	POP-START
		name		SouthernWV
		samples		WVA_01 d WVA_02 d WVA_06 d WVA_20 d		
	POP-END

	POP-START
		name		Bluegrass
		samples		EKY_02 d EKY_03 d EKY_06 d EKY_04 d EKY_01 d
	POP-END

	POP-START
		name 		NorthernWV
		samples		WVA_03 d WVA_07 d WVA_08 d WVA_11 d 	
	POP-END

	
CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			Appalachia1
		children		Smokie		SouthernWV
		tau-initial		0.00001
	POP-END
	
	POP-START
		name			BlueAppalachia1
		children		Bluegrass		Appalachia1
		tau-initial		0.00001
	POP-END
	
	POP-START
		name			Root
		children		NorthernWV		BlueAppalachia1
		tau-initial		0.00001
	POP-END

ANCESTRAL-POPS-END

MIG-BANDS-START	

	BAND-START		
       source  NorthernWV
       target  Bluegrass
	BAND-END
	
	BAND-START		
       source  NorthernWV
       target  Smokie
	BAND-END
	
	BAND-START		
       source  NorthernWV
       target  SouthernWV
	BAND-END
	
	BAND-START		
       source  Bluegrass
       target  Smokie
	BAND-END
	
	BAND-START		
       source  Bluegrass
       target  NorthernWV
	BAND-END
	
	BAND-START		
       source  Bluegrass
       target  SouthernWV
	BAND-END
	
	BAND-START		
       source  Smokie
       target  Bluegrass
	BAND-END
	
	BAND-START		
       source  Smokie
       target  NorthernWV
	BAND-END
	
	BAND-START		
       source  Smokie
       target  SouthernWV
    BAND-END
	
	BAND-START		
       source  SouthernWV
       target  Bluegrass
    BAND-END

	BAND-START		
       source  SouthernWV
       target  Smokie
    BAND-END
	
	BAND-START		
       source  SouthernWV
       target  NorthernWV
    BAND-END

	BAND-START
	source Appalachia1
	target Bluegrass
    BAND-END

	BAND-START
	source BlueAppalachia1
	target NorthernWV
    BAND-END

MIG-BANDS-END
