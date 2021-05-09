from sage.libs.pari.convert_sage import gen_to_sage
from sage.rings.number_field.splitting_field import SplittingData

def LLLbasiscore(L):
	L.inject_variables()
	Com = ComplexField(10000)
	Rf = RealField(10000)

	bnf = L.pari_bnf(proof = False, units=True)
	units = bnf[7][4]
	embeddings = L.embeddings(Com)
	embeds = len(embeddings)
	funds = len(units)
	print("bnf done")
	
	basis = []
	for u in units:
		basis.append([log(abs(e(u))) for e in embeddings])
		
	basis0 = []
	for i in range(funds):
		basis0.append([])
		for j in range(embeds):
			basis0[i].append(basis[i][j])
	basis0.append([])
	for j in range(embeds):
		basis0[funds].append(1)
	N = Matrix(Rf, funds+1, embeds, basis0)
	determ = Rf(N.determinant())/Rf(sqrt(embeds))
	
	large = 10**(100)
	for i in range(funds):
		for j in range(embeds):
			basis[i][j] *= large
			basis[i][j] = basis[i][j].integer_part();
	M = Matrix(ZZ, funds, embeds, basis)
	prod0 = Rf(1)
	for i in range(funds):
		mag1 = 0
		for j in range(embeds):
			mag1 += Rf(N[i][j])^2
		prod0 *= mag1
	prod0 = Rf(sqrt(prod0))
	
	LLLtemp = M.LLL(eta=0.51)
	prodL = Rf(1)
	mag = 0
	basisL = []
	shortL=0
	for i in range(funds):
		basisL.append([])
		mag = 0
		for j in range(embeds):
			basisL[i].append(Rf(LLLtemp[i][j])/large)
			mag += (Rf(LLLtemp[i][j])/large)^2
		mag = sqrt(mag)
		#print(RR(mag))
		prodL *= mag
		if i == 0:
			shortL = mag
	LLL = Matrix(RR, funds,embeds, basisL)
	print("LLL done")
	
	BKZtemp = M.BKZ()
	prodB = Rf(1)
	mag = 0
	basisB = []
	shortB = 0
	for i in range(funds):
		basisB.append([])
		mag = 0
		for j in range(embeds):
			basisB[i].append(RR(BKZtemp[i][j])/large)
			mag += (Rf(BKZtemp[i][j])/large)^2
		mag = sqrt(mag)
		#print(RR(mag))
		prodB *= mag
		if i == 0:
			shortB = mag
	BKZ = Matrix(RR, funds,embeds, basisB)
	print("bkz done")
	
	return [embeds, N.change_ring(RR), LLL, BKZ,[prod0,prodL,prodB], abs(determ), [shortL,shortB]]

def LLLpoly(p):
	try:
		K.<a> = NumberField(p)
		print(str(p) + " in progress")
		return [p,LLLbasiscore(K)]
	except:
		print("Something is wrong with the polynomial:")
		print(p)
		return []
	
def hilbClassLLL(n):
	gpoly = gp.quadhilbert(n)
	firstpoly = firstDeg(gp.polcompositum(gpoly,gpoly),24)
	if firstpoly == 0:
		print(str(n) + " messed up somehow")
	genpoly = gen_to_sage(pari(firstpoly),{'x':x})
	return LLLpoly(genpoly)

def cyclotomicLLL(n):
	gpoly = gp.nfsubfields(gp.polcyclo(n),24)[1][1]
	gpoly = gp.polredbest(gpoly)
	genpoly = gen_to_sage(pari(gpoly),{'x':x})
	return LLLpoly(genpoly)

def cyclicPairLLL(p,q,n,m):
	gpolyp = gp.nfsubfields(gp.polcyclo(p),n)[1][1]
	gpolyp = gp.polredbest(gpolyp)
	gpolyq = gp.nfsubfields(gp.polcyclo(q),m)[1][1]
	gpolyq = gp.polredbest(gpolyq)
	print(gpolyp)
	print(gpolyq)
	firstpoly = firstDeg(gp.polcompositum(gpolyp,gpolyq),24)
	if firstpoly == 0:
		print(str(p) + "," + str(q) + " messed up somehow")
	else:
		genpoly = gen_to_sage(pari(firstpoly),{'x':x})
		return LLLpoly(genpoly)

def splitLLL(poly):
	L.<a> = poly.splitting_field()
	print(str(poly) + " in progress")
	return [poly,LLLbasiscore(L)]

def checkQuartic(poly):
	if ((gp.polisirreducible(poly) == 1) and (gp.polgalois(poly)[1]==24)):
		L.<a> = NumberField(poly)
		return (len(L.embeddings(RR))==4)
	return 0
		

def searchQuartic(bound):
	goodQuartics = []
	for i in range(-1*bound,bound+1):
		for j in range(-1*bound,bound+1):
			for k in range(-1*bound,bound+1):
				for l in range(-1*bound,bound+1):
					p = x^4 + i*x^3 + j*x^2 + k*x+l
					if checkQuartic(p):
						print(p)
						goodQuartics.append(p)
	return goodQuartics
	
def parseCore(coreData):
	#print("LLL Reduction Factor: \t\t" + str(RR(coreData[4][0]/coreData[4][1])))
	print("LLL Hermite Defect: \t\t" + str(RR(coreData[6][0]/(coreData[5]**(1/coreData[0])))))
	print("LLL Orthogonality Defect: \t" + str(RR(coreData[4][1]/coreData[5])))
	#print("BKZ Reduction Factor: \t\t" + str(RR(coreData[4][0]/coreData[4][2])))
	print("BKZ Hermite Defect: \t\t" + str(RR(coreData[6][1]/(coreData[5]**(1/coreData[0])))))
	print("BKZ Orthogonality Defect: \t" + str(RR(coreData[4][2]/coreData[5])))
	return [RR(coreData[6][0]/(coreData[5]**(1/coreData[0]))), RR(coreData[4][1]/coreData[5]), RR(coreData[6][1]/(coreData[5]**(1/coreData[0]))), RR(coreData[4][2]/coreData[5])]
	
def hilbClassParse(n):
	gpoly = gp.quadhilbert(n)
	firstpoly = firstDeg(gp.polcompositum(gpoly,gpoly),24)
	if firstpoly == 0:
		print(str(n) + " messed up somehow")
	genpoly = gen_to_sage(pari(firstpoly),{'x':x})
	K.<a> = NumberField(genpoly)
	print(K)
	parseCore(LLLbasiscore(K))

def parseCoreList(dataList):
	relevantData = []
	relevantDataTranspose = [[],[],[],[]]
	for data in dataList:
		temp = parseCore(data[1])
		relevantData.append(temp)
		for j in range(4):
			relevantDataTranspose[j].append(temp[j])
		print("")
	return [relevantData,relevantDataTranspose]
	
def firstDeg(gpPolys, n):
	for p in gpPolys:
		if gp.poldegree(p)==n:
			return p
	return 0;
	
def dataDump(quartics):
	C6nums = C6List()
	C12nums = C12List()
	C24nums = C24List()
	C8nums = C24p1List()
	C3nums = C24p2List()
	S4quartics = quartics
	C6data = []
	C12data = []
	C24data = []
	S4data = []

	for n in C6nums:
		C6data.append(hilbClassLLL(n))
	for n in C12nums:
		C12data.append(hilbClassLLL(n))
	for n in C24nums:
		C24data.append(cyclotomicLLL(n))
	for i in range(len(C8nums)):
		C24data.append(cyclicPairLLL(C8nums[i],C3nums[i],8,3))
	for poly in quartics:
		S4data.append(splitLLL(poly))
		
	return [[C6data,C6nums],[C12data,C12nums],[C24data,C24nums,C8nums,C3nums],[S4data,quartics]]

def dumpNonCyc(quartics):
	C6nums = C6List()
	C12nums = C12List()
	S4quartics = quartics
	C6data = []
	C12data = []
	S4data = []
	
	for n in C6nums:
		C6data.append(hilbClassLLL(n))
	for n in C12nums:
		C12data.append(hilbClassLLL(n))
	for poly in quartics:
		S4data.append(splitLLL(poly))

	return [[C6data,C6nums],[C12data,C12nums],[S4data,quartics]]

def dumpCyclic():
	C24nums = C24List()
	C8nums = C24p1List()
	C3nums = C24p2List()
	C24data = []
	for n in C24nums:
		C24data.append(cyclotomicLLL(n))
	for i in range(len(C8nums)):
		C24data.append(cyclicPairLLL(C8nums[i],C3nums[i],8,3))
	return C24data

def jankyCyclic(n8,n3):
	C24nums = C24List()
	C8nums = C24p1List()
	C3nums = C24p2List()
	C24data = []
	for n in C24nums:
		C24data.append(cyclotomicLLL(n))
	for i in range(n8):
		for j in range(n3):
			C24data.append(cyclicPairLLL(C8nums[i],C3nums[j],8,3))
	return C24data
	
#ugly stuff/useful variable lists
		
def polyList():
	return [x^10 - 15*x^8 + x^7 + 49*x^6 - 50*x^4 + 16*x^2 - 1,
	x^10 - 15*x^8 + x^7 + 50*x^6 - 42*x^4 + 12*x^2 - 1,
	x^10 - 15*x^8 + x^7 + 50*x^6 - 43*x^4 + 12*x^2 - 1,
	x^10 - 15*x^8 + x^7 + 50*x^6 - 44*x^4 + 12*x^2 - 1,
	x^10 - 15*x^8 + x^7 + 50*x^6 - 44*x^4 + 13*x^2 - 1,
	x^10 - 15*x^8 + x^7 + 50*x^6 - 45*x^4 + 13*x^2 - 1,
	x^10 - 15*x^8 + x^7 + 50*x^6 - 46*x^4 + 13*x^2 - 1,
	x^10 - 15*x^8 + x^7 + 50*x^6 - 46*x^4 + 14*x^2 - 1,
	x^10 - 16*x^8 + x^7 + 45*x^6 - 37*x^4 + 11*x^2 - 1,
	x^10 - 16*x^8 + x^7 + 45*x^6 - 38*x^4 + 11*x^2 - 1,
	x^10 - 16*x^8 + x^7 + 45*x^6 - 39*x^4 + 12*x^2 - 1,
	x^10 - 12*x^8 + x^7 + 38*x^6 - 39*x^4 + 13*x^2 - 1,
	x^10 - 12*x^8 + x^7 + 38*x^6 - 40*x^4 + 14*x^2 - 1,
	x^10 - 12*x^8 + x^7 + 39*x^6 - 35*x^4 + 11*x^2 - 1,
	x^10 - 12*x^8 + x^7 + 39*x^6 - 36*x^4 + 11*x^2 - 1]

def C12List():
	return [4345,
	5629,
	6088,
	6401,
	6856,
	7745,
	7881,
	8545,
	9676,
	11656,
	13549,
	13801,
	14424,
	15897,
	16409,
	16609,
	16645,
	16684,
	16913,
	17113,
	17884,
	17905,
	18409,
	19113,
	19756,
	20168,
	21064,
	21324,
	21641,
	22476,
	22633,
	22873,
	23377,
	23953,
	24145,
	24433,
	27224,
	27289,
	27445,
	27713,
	27833,
	28024,
	28041,
	28381,
	29245,
	29905,
	30561,
	31001,
	31432,
	32097]

def C6List():
	return [2920,
	4360,
	4641,
	5980,
	6396,
	7084,
	7224,
	7665,
	8220,
	8556,
	8680,
	8745,
	8905,
	9805,
	10353,
	10540,
	10812,
	10865,
	11020,
	11505,
	12441,
	12765,
	14385,
	14824,
	15820,
	15873,
	16044,
	16385,
	16860,
	17420,
	18492,
	18824,
	19020,
	19045,
	19880,
	19885,
	20865,
	22044,
	22380,
	23720,
	23816,
	24232,
	24645,
	24648,
	25185,
	25228,
	25705,
	26376,
	26680,
	26845]

def C24List():
	return [97,
	193,
	241,
	337,
	433,
	577,
	673,
	769,
	1009,
	1153,
	1201,
	1249,
	1297,
	1489,
	1777,
	1873,
	2017]

def C24p1List():
	return [113,
	257,
	353,
	401,
	449,
	593,
	641,
	881,
	929,
	977,
	1217,
	1361,
	1409,
	1553,
	1601,
	1697,
	1889,
	2081,
	2129,
	2273,
	2417,
	2609,
	2657,
	2753,
	2801,
	2897,
	3041,
	3089,
	3137,
	3329,
	3617]

def C24p2List():
	return [13,
	19,
	31,
	37,
	43,
	61,
	67,
	79,
	103,
	109,
	127,
	139,
	151,
	157,
	163,
	181,
	199,
	211,
	223,
	229,
	271,
	277,
	283,
	307,
	331,
	349,
	367,
	373,
	379,
	397,
	421]
