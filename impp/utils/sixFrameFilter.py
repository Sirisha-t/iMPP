import math
import re
import sys
import string
import time


def CodonTranslation():

	global  gl_Ter
	global  gl_Phe, gl_Leu, gl_Ser, gl_Tyr, gl_Cys, \
		gl_Trp, gl_Pro, gl_His, gl_Gln, gl_Arg, \
		gl_Ile, gl_Met, gl_Thr, gl_Asn, gl_Lys, \
		gl_Val, gl_Ala, gl_Asp, gl_Glu, gl_Gly
	global  gl_Any

	global  gl_AmAcLib

	gl_Ter = [ "TAA", "TAG", "TGA" ]

	gl_Phe = [ "TTT", "TTC" ]
	gl_Leu = [ "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "CTN" ]
	gl_Ser = [ "TCT", "TCC", "TCA", "TCG", "TCN", "AGT", "AGC" ]
	gl_Tyr = [ "TAT", "TAC" ]
	gl_Cys = [ "TGT", "TGC" ]

	gl_Trp = [ "TGG" ]
	gl_Pro = [ "CCT", "CCC", "CCA", "CCG", "CCN" ]
	gl_His = [ "CAT", "CAC" ]
	gl_Gln = [ "CAA", "CAG" ]
	gl_Arg = [ "CGT", "CGC", "CGA", "CGG", "CGN", "AGA", "AGG" ]

	gl_Ile = [ "ATT", "ATC", "ATA" ]
	gl_Met = [ "ATG" ]
	gl_Thr = [ "ACT", "ACC", "ACA", "ACG", "ACN" ]
	gl_Asn = [ "AAT", "AAC" ]
	gl_Lys = [ "AAA", "AAG" ]

	gl_Val = [ "GTT", "GTC", "GTA", "GTG", "GTN" ]
	gl_Ala = [ "GCT", "GCC", "GCA", "GCG", "GCN" ]
	gl_Asp = [ "GAT", "GAC" ]
	gl_Glu = [ "GAA", "GAG" ]
	gl_Gly = [ "GGT", "GGC", "GGA", "GGG", "GGN" ]

	gl_Any = [ "AAN", "ATN", "AGN", \
		   "TAN", "TTN", "TGN", \
		   "GAN", \
		   "CAN", \
		   "ANA", "ANT", "ANG", "ANC", \
		   "TNA", "TNT", "TNG", "TNC", \
		   "GNA", "GNT", "GNG", "GNC", \
		   "CNA", "CNT", "CNG", "CNC", \
		   "NAA", "NAT", "NAG", "NAC", \
		   "NTA", "NTT", "NTG", "NTC", \
		   "NGA", "NGT", "NGG", "NGC", \
		   "NCA", "NCT", "NCG", "NCC", \
		   "NNA", "NNT", "NNG", "NNC", \
		   "ANN", "TNN", "GNN", "CNN", \
		   "NAN", "NTN", "NGN", "NCN", \
		   "NNN" ]
	######################################
	gl_AmAcLib =  [ gl_Ter, \
			gl_Phe, gl_Leu, gl_Ser, gl_Tyr, gl_Cys, \
			gl_Trp, gl_Pro, gl_His, gl_Gln, gl_Arg, \
			gl_Ile, gl_Met, gl_Thr, gl_Asn, gl_Lys, \
			gl_Val, gl_Ala, gl_Asp, gl_Glu, gl_Gly, \
			gl_Any ]
	######################################

def NSeq_Processor(have_seqs):
	non_atgc_list = [ 'B', 'D', 'E', 'F', 'H', 'I', 'J', 'K', 'L', 'M', \
			'O', 'P', 'Q', 'R', 'S', 'U', 'V', 'W', 'X', 'Y', 'Z' ]

	have_seqs = string.upper(have_seqs)
	for dummy_letter in non_atgc_list:
		have_seqs = re.sub(dummy_letter, "N", have_seqs)

	return have_seqs

def Seqs_Translator(proper_id, have_seqs, trans_fr, gen_code):
	abc_list =     ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', \
			'J','K','L','M','N','O','P','Q','R','S','T', \
			'U','V','W','X','Y','Z']
	bad_list =     [     'B',      'D', 'E', 'F',      'H', 'I', \
			'J','K','L','M',    'O','P','Q','R','S',     \
			'U','V','W','X','Y','Z'] # REMOVED A T G C N

        global  gl_Ter
	global  gl_Phe, gl_Leu, gl_Ser, gl_Tyr, gl_Cys, \
		gl_Trp, gl_Pro, gl_His, gl_Gln, gl_Arg, \
		gl_Ile, gl_Met, gl_Thr, gl_Asn, gl_Lys, \
		gl_Val, gl_Ala, gl_Asp, gl_Glu, gl_Gly
	global  gl_Any

	global  gl_AmAcLib

	global  gl_longest_orf

        t = have_seqs			# DNA sequence
	t = string.upper(t)		# all uppercase
	t = re.sub("-", "", t)		# remove all dashes
	t = list(t)			# string -> list of chars
	gl_longest_orf = ""
        mooba = 0	# COUNTER OF ALL LETTERS
	booba = 0	# COUNTER OF BAD LETTERS
	for item in t:
		if item in bad_list:
			t[mooba] = "N"
			booba = booba + 1
		mooba = mooba + 1
        #########################################################
 	t_r = t[:]			# will be reverse string
	t_r.reverse()			# inplace reverse the list
	t   = string.join(t,   '')	# list of strings -> forward string
	t_r = string.join(t_r, '')	# list of strings -> reverse string
	#########################################
	t_r = re.sub("A", "t", t_r)	# A -> t
	t_r = re.sub("T", "a", t_r)	# T -> a
	t_r = re.sub("G", "c", t_r)	# G -> c
	t_r = re.sub("C", "g", t_r)	# C -> g
	#########################################
	t_r = string.upper(t_r)	# back to uppercase
	l = len(t)
  	 ### TRANSLATION SIX FRAMES ###
	if trans_fr == 6:
		k1 = 0
		k2 = 1
		k3 = 2
		###############
		stop_count1 = 0
		any_count1  = 0
		longest_fr1 = 0
		longest_orf1 = ""
		total_size1 = 0
		frm_number1 = 0
		my_list1 = []
		###############
		stop_count2 = 0
		any_count2  = 0
		longest_fr2 = 0
		longest_orf2 = ""
		total_size2 = 0
		frm_number2 = 0
		my_list2 = []
		###############
		stop_count3 = 0
		any_count3  = 0
		longest_fr3 = 0
		longest_orf3 = ""
		total_size3 = 0
		frm_number3 = 0
		my_list3 = []
		###############
		stop_count4 = 0
		any_count4  = 0
		longest_fr4 = 0
		longest_orf4 = ""
		total_size4 = 0
		frm_number4 = 0
		my_list4 = []
		###############
		stop_count5 = 0
		any_count5  = 0
		longest_fr5 = 0
		longest_orf5 = ""
		total_size5 = 0
		frm_number5 = 0
		my_list5 = []
		###############
		stop_count6 = 0
		any_count6  = 0
		longest_fr6 = 0
		longest_orf6 = ""
		total_size6 = 0
		frm_number6 = 0
		my_list6 = []
		###############
		while k1 < l:
			m1 = k1 + 3
			m2 = k2 + 3
			m3 = k3 + 3
			p1 = t[k1:m1]	# triplet frame 1
			p2 = t[k2:m2]	# triplet frame 2
			p3 = t[k3:m3]	# triplet frame 3
			p4 = t_r[k1:m1]	# triplet frame 4
			p5 = t_r[k2:m2]	# triplet frame 5
			p6 = t_r[k3:m3]	# triplet frame 6
			#################################
			if p1 in gl_Any:
				my_list1.append("X")
				any_count1  =  any_count1 + 1
			if p1 in gl_Ter:
				my_list1.append("*")
				stop_count1 = stop_count1 + 1
			if p1 in gl_Phe:
				my_list1.append("F")
			if p1 in gl_Leu:
				my_list1.append("L")
			if p1 in gl_Ser:
				my_list1.append("S")
			if p1 in gl_Tyr:
				my_list1.append("Y")
			if p1 in gl_Cys:
				my_list1.append("C")
			if p1 in gl_Trp:
				my_list1.append("W")
			if p1 in gl_Pro:
				my_list1.append("P")
			if p1 in gl_His:
				my_list1.append("H")
			if p1 in gl_Gln:
				my_list1.append("Q")
			if p1 in gl_Arg:
				my_list1.append("R")
			if p1 in gl_Ile:
				my_list1.append("I")
			if p1 in gl_Met:
				my_list1.append("M")
			if p1 in gl_Thr:
				my_list1.append("T")
			if p1 in gl_Asn:
				my_list1.append("N")
			if p1 in gl_Lys:
				my_list1.append("K")
			if p1 in gl_Val:
				my_list1.append("V")
			if p1 in gl_Ala:
				my_list1.append("A")
			if p1 in gl_Asp:
				my_list1.append("D")
			if p1 in gl_Glu:
				my_list1.append("E")
			if p1 in gl_Gly:
				my_list1.append("G")
			##################################
			if p2 in gl_Any:
				my_list2.append("X")
				any_count2  =  any_count2 + 1
			if p2 in gl_Ter:
				my_list2.append("*")
				stop_count2 = stop_count2 + 1
			if p2 in gl_Phe:
				my_list2.append("F")
			if p2 in gl_Leu:
				my_list2.append("L")
			if p2 in gl_Ser:
				my_list2.append("S")
			if p2 in gl_Tyr:
				my_list2.append("Y")
			if p2 in gl_Cys:
				my_list2.append("C")
			if p2 in gl_Trp:
				my_list2.append("W")
			if p2 in gl_Pro:
				my_list2.append("P")
			if p2 in gl_His:
				my_list2.append("H")
			if p2 in gl_Gln:
				my_list2.append("Q")
			if p2 in gl_Arg:
				my_list2.append("R")
			if p2 in gl_Ile:
				my_list2.append("I")
			if p2 in gl_Met:
				my_list2.append("M")
			if p2 in gl_Thr:
				my_list2.append("T")
			if p2 in gl_Asn:
				my_list2.append("N")
			if p2 in gl_Lys:
				my_list2.append("K")
			if p2 in gl_Val:
				my_list2.append("V")
			if p2 in gl_Ala:
				my_list2.append("A")
			if p2 in gl_Asp:
				my_list2.append("D")
			if p2 in gl_Glu:
				my_list2.append("E")
			if p2 in gl_Gly:
				my_list2.append("G")
			##################################
			if p3 in gl_Any:
				my_list3.append("X")
				any_count3  =  any_count3 + 1
			if p3 in gl_Ter:
				my_list3.append("*")
				stop_count3 = stop_count3 + 1
			if p3 in gl_Phe:
				my_list3.append("F")
			if p3 in gl_Leu:
				my_list3.append("L")
			if p3 in gl_Ser:
				my_list3.append("S")
			if p3 in gl_Tyr:
				my_list3.append("Y")
			if p3 in gl_Cys:
				my_list3.append("C")
			if p3 in gl_Trp:
				my_list3.append("W")
			if p3 in gl_Pro:
				my_list3.append("P")
			if p3 in gl_His:
				my_list3.append("H")
			if p3 in gl_Gln:
				my_list3.append("Q")
			if p3 in gl_Arg:
				my_list3.append("R")
			if p3 in gl_Ile:
				my_list3.append("I")
			if p3 in gl_Met:
				my_list3.append("M")
			if p3 in gl_Thr:
				my_list3.append("T")
			if p3 in gl_Asn:
				my_list3.append("N")
			if p3 in gl_Lys:
				my_list3.append("K")
			if p3 in gl_Val:
				my_list3.append("V")
			if p3 in gl_Ala:
				my_list3.append("A")
			if p3 in gl_Asp:
				my_list3.append("D")
			if p3 in gl_Glu:
				my_list3.append("E")
			if p3 in gl_Gly:
				my_list3.append("G")
			##################################
			if p4 in gl_Any:
				my_list4.append("X")
				any_count4  =  any_count4 + 1
			if p4 in gl_Ter:
				my_list4.append("*")
				stop_count4 = stop_count4 + 1
			if p4 in gl_Phe:
				my_list4.append("F")
			if p4 in gl_Leu:
				my_list4.append("L")
			if p4 in gl_Ser:
				my_list4.append("S")
			if p4 in gl_Tyr:
				my_list4.append("Y")
			if p4 in gl_Cys:
				my_list4.append("C")
			if p4 in gl_Trp:
				my_list4.append("W")
			if p4 in gl_Pro:
				my_list4.append("P")
			if p4 in gl_His:
				my_list4.append("H")
			if p4 in gl_Gln:
				my_list4.append("Q")
			if p4 in gl_Arg:
				my_list4.append("R")
			if p4 in gl_Ile:
				my_list4.append("I")
			if p4 in gl_Met:
				my_list4.append("M")
			if p4 in gl_Thr:
				my_list4.append("T")
			if p4 in gl_Asn:
				my_list4.append("N")
			if p4 in gl_Lys:
				my_list4.append("K")
			if p4 in gl_Val:
				my_list4.append("V")
			if p4 in gl_Ala:
				my_list4.append("A")
			if p4 in gl_Asp:
				my_list4.append("D")
			if p4 in gl_Glu:
				my_list4.append("E")
			if p4 in gl_Gly:
				my_list4.append("G")
			##################################
			if p5 in gl_Any:
				my_list5.append("X")
				any_count5  =  any_count5 + 1
			if p5 in gl_Ter:
				my_list5.append("*")
				stop_count5 = stop_count5 + 1
			if p5 in gl_Phe:
				my_list5.append("F")
			if p5 in gl_Leu:
				my_list5.append("L")
			if p5 in gl_Ser:
				my_list5.append("S")
			if p5 in gl_Tyr:
				my_list5.append("Y")
			if p5 in gl_Cys:
				my_list5.append("C")
			if p5 in gl_Trp:
				my_list5.append("W")
			if p5 in gl_Pro:
				my_list5.append("P")
			if p5 in gl_His:
				my_list5.append("H")
			if p5 in gl_Gln:
				my_list5.append("Q")
			if p5 in gl_Arg:
				my_list5.append("R")
			if p5 in gl_Ile:
				my_list5.append("I")
			if p5 in gl_Met:
				my_list5.append("M")
			if p5 in gl_Thr:
				my_list5.append("T")
			if p5 in gl_Asn:
				my_list5.append("N")
			if p5 in gl_Lys:
				my_list5.append("K")
			if p5 in gl_Val:
				my_list5.append("V")
			if p5 in gl_Ala:
				my_list5.append("A")
			if p5 in gl_Asp:
				my_list5.append("D")
			if p5 in gl_Glu:
				my_list5.append("E")
			if p5 in gl_Gly:
				my_list5.append("G")
			##################################
			if p6 in gl_Any:
				my_list6.append("X")
				any_count6  =  any_count6 + 1
			if p6 in gl_Ter:
				my_list6.append("*")
				stop_count6 = stop_count6 + 1
			if p6 in gl_Phe:
				my_list6.append("F")
			if p6 in gl_Leu:
				my_list6.append("L")
			if p6 in gl_Ser:
				my_list6.append("S")
			if p6 in gl_Tyr:
				my_list6.append("Y")
			if p6 in gl_Cys:
				my_list6.append("C")
			if p6 in gl_Trp:
				my_list6.append("W")
			if p6 in gl_Pro:
				my_list6.append("P")
			if p6 in gl_His:
				my_list6.append("H")
			if p6 in gl_Gln:
				my_list6.append("Q")
			if p6 in gl_Arg:
				my_list6.append("R")
			if p6 in gl_Ile:
				my_list6.append("I")
			if p6 in gl_Met:
				my_list6.append("M")
			if p6 in gl_Thr:
				my_list6.append("T")
			if p6 in gl_Asn:
				my_list6.append("N")
			if p6 in gl_Lys:
				my_list6.append("K")
			if p6 in gl_Val:
				my_list6.append("V")
			if p6 in gl_Ala:
				my_list6.append("A")
			if p6 in gl_Asp:
				my_list6.append("D")
			if p6 in gl_Glu:
				my_list6.append("E")
			if p6 in gl_Gly:
				my_list6.append("G")
			##################################
			k1 = k1 + 3
			k2 = k2 + 3	# next triplet frame2
			k3 = k3 + 3	# next triplet frame3
			###########
                my_string1 = "".join(my_list1)
		my_string2 = "".join(my_list2)
		my_string3 = "".join(my_list3)
		my_string4 = "".join(my_list4)
		my_string5 = "".join(my_list5)
		my_string6 = "".join(my_list6)
		#############################
		total_size1 = len(my_string1)
		total_size2 = len(my_string2)
		total_size3 = len(my_string3)
		total_size4 = len(my_string4)
		total_size5 = len(my_string5)
		total_size6 = len(my_string6)
		#############################
		my_subset1 = my_string1.split('*')
		my_subset2 = my_string2.split('*')
		my_subset3 = my_string3.split('*')
		my_subset4 = my_string4.split('*')
		my_subset5 = my_string5.split('*')
		my_subset6 = my_string6.split('*')

        	## 1
		for fragment1 in my_subset1:
			fragment_len1 = len(fragment1)
			# if fragment_len1 >= longest_fr1:
			if fragment_len1 > longest_fr1:
				longest_fr1 = fragment_len1
				longest_orf1 = fragment1
			if fragment_len1 != 0:
				frm_number1 = frm_number1 + 1
		## 2
		for fragment2 in my_subset2:
			fragment_len2 = len(fragment2)
			# if fragment_len2 >= longest_fr2:
			if fragment_len2 > longest_fr2:
				longest_fr2 = fragment_len2
				longest_orf2 = fragment2
			if fragment_len2 != 0:
				frm_number2 = frm_number2 + 1
		## 3
		for fragment3 in my_subset3:
			fragment_len3 = len(fragment3)
			# if fragment_len3 >= longest_fr3:
			if fragment_len3 > longest_fr3:
				longest_fr3 = fragment_len3
				longest_orf3 = fragment3
			if fragment_len3 != 0:
				frm_number3 = frm_number3 + 1
		## 4
		for fragment4 in my_subset4:
			fragment_len4 = len(fragment4)
			# if fragment_len4 >= longest_fr4:
			if fragment_len4 > longest_fr4:
				longest_fr4 = fragment_len4
				longest_orf4 = fragment4
			if fragment_len4 != 0:
				frm_number4 = frm_number4 + 1
		## 5
		for fragment5 in my_subset5:
			fragment_len5 = len(fragment5)
			# if fragment_len5 >= longest_fr5:
			if fragment_len5 > longest_fr5:
				longest_fr5 = fragment_len5
				longest_orf5 = fragment5
			if fragment_len5 != 0:
				frm_number5 = frm_number5 + 1
		## 6
		for fragment6 in my_subset6:
			fragment_len6 = len(fragment6)
			# if fragment_len6 >= longest_fr6:
			if fragment_len6 > longest_fr6:
				longest_fr6 = fragment_len6
				longest_orf6 = fragment6
			if fragment_len6 != 0:
				frm_number6 = frm_number6 + 1
        	############################
        	full_frame1 = "UNDEFINED"
		full_frame2 = "UNDEFINED"
		full_frame3 = "UNDEFINED"
		full_frame4 = "UNDEFINED"
		full_frame5 = "UNDEFINED"
		full_frame6 = "UNDEFINED"
		#############################
        	## 1
		if longest_fr1 >= total_size1 - 1:
			full_frame1 = "SINGLE_ORF"
		if longest_fr1 < total_size1 - 1:
			full_frame1 = "MULTIPLE_ORFs"
		## 2
		if longest_fr2 >= total_size2 - 1:
			full_frame2 = "SINGLE_ORF"
		if longest_fr2 < total_size2 - 1:
			full_frame2 = "MULTIPLE_ORFs"
		## 3
		if longest_fr3 >= total_size3 - 1:
			full_frame3 = "SINGLE_ORF"
		if longest_fr3 < total_size3 - 1:
			full_frame3 = "MULTIPLE_ORFs"
		## 4
		if longest_fr4 >= total_size4 - 1:
			full_frame4 = "SINGLE_ORF"
		if longest_fr4 < total_size4 - 1:
			full_frame4 = "MULTIPLE_ORFs"
		## 5
		if longest_fr5 >= total_size5 - 1:
			full_frame5 = "SINGLE_ORF"
		if longest_fr5 < total_size5 - 1:
			full_frame5 = "MULTIPLE_ORFs"
		## 6
		if longest_fr6 >= total_size6 - 1:
			full_frame6 = "SINGLE_ORF"
		if longest_fr6 < total_size6 - 1:
			full_frame6 = "MULTIPLE_ORFs"
		#############################
        	## LONGEST ##
		llen1 = len(longest_orf1)
		llen2 = len(longest_orf2)
		llen3 = len(longest_orf3)
		llen4 = len(longest_orf4)
		llen5 = len(longest_orf5)
		llen6 = len(longest_orf6)
		lfrm  = " "
		lfrm1 = ""
		lfrm2 = ""
		lfrm3 = ""
		lfrm4 = ""
		lfrm5 = ""
		lfrm6 = ""
		dupl_status = " UNIQ_LONG"
		## 1
		if llen1 >= llen2 and llen1 >= llen3 and llen1 >= llen4 and llen1 >= llen5 and llen1 >= llen6:
			if llen1 == llen2 or llen1 == llen3 or llen1 == llen4 or llen1 == llen5 or llen1 == llen6:
				dupl_status = " DUPL_LONG"
			lfr = "_fr1"
			lfrm1 = "fr_1"
			longest_orf = longest_orf1
			stop_count = stop_count1
			any_count = any_count1
			longest_fr = longest_fr1
			total_size = total_size1
			full_frame = full_frame1
			frm_number = frm_number1
			longest_len = llen1
			#gl_longest_orf.write(">" + proper_id + lfr + dupl_status + \
			#		" STOP:" + `stop_count` + " ANY(X):" + `any_count` + \
			#		" LONGEST:" + `longest_fr` + " TOTAL_LENGTH:" + `total_size` + \
			#		" " + full_frame + ":" + `frm_number` + '\n')
			#print longest_orf
			gl_longest_orf += ">" + proper_id + lfr + '\n' + longest_orf + '\n'
			#gl_longest_orf.write(longest_orf + '\n')
		## 2
		if llen2 >= llen1 and llen2 >= llen3 and llen2 >= llen4 and llen2 >= llen5 and llen2 >= llen6:
			if llen2 == llen1 or llen2 == llen3 or llen2 == llen4 or llen2 == llen5 or llen2 == llen6:
				dupl_status = " DUPL_LONG"
			lfr = "_fr2"
			lfrm2 = "fr_2"
			longest_orf = longest_orf2
			stop_count = stop_count2
			any_count = any_count2
			longest_fr = longest_fr2
			total_size = total_size2
			full_frame = full_frame2
			frm_number = frm_number2
			longest_len = llen2
			'''
			gl_longest_orf.write(">" + proper_id + lfr + dupl_status + \
					" STOP:" + `stop_count` + " ANY(X):" + `any_count` + \
					" LONGEST:" + `longest_fr` + " TOTAL_LENGTH:" + `total_size` + \
					" " + full_frame + ":" + `frm_number` + '\n')
			#print longest_orf
			gl_longest_orf.write(longest_orf + '\n')
			'''
			gl_longest_orf += ">" + proper_id + lfr + '\n' + longest_orf + '\n'
		## 3
		if llen3 >= llen1 and llen3 >= llen2 and llen3 >= llen4 and llen3 >= llen5 and llen3 >= llen6:
			if llen3 == llen1 or llen3 == llen2 or llen3 == llen4 or llen3 == llen5 or llen3 == llen6:
				dupl_status = " DUPL_LONG"
			lfr = "_fr3"
			lfrm3 = "fr_3"
			longest_orf = longest_orf3
			stop_count = stop_count3
			any_count = any_count3
			longest_fr = longest_fr3
			total_size = total_size3
			full_frame = full_frame3
			frm_number = frm_number3
			longest_len = llen3
			'''
			gl_longest_orf.write(">" + proper_id + lfr + dupl_status + \
					" STOP:" + `stop_count` + " ANY(X):" + `any_count` + \
					" LONGEST:" + `longest_fr` + " TOTAL_LENGTH:" + `total_size` + \
					" " + full_frame + ":" + `frm_number` + '\n')
			#print longest_orf
			gl_longest_orf.write(longest_orf + '\n')
			'''
			gl_longest_orf += ">" + proper_id + lfr + '\n' + longest_orf + '\n'
		## 4	
		if llen4 >= llen1 and llen4 >= llen2 and llen4 >= llen3 and llen4 >= llen5 and llen4 >= llen6:
			if llen4 == llen1 or llen4 == llen2 or llen4 == llen3 or llen4 == llen5 or llen4 == llen6:
				dupl_status = " DUPL_LONG"
			lfr = "_fr4"
			lfrm4 = "fr_4"
			longest_orf = longest_orf4
			stop_count = stop_count4
			any_count = any_count4
			longest_fr = longest_fr4
			total_size = total_size4
			full_frame = full_frame4
			frm_number = frm_number4
			longest_len = llen4
			'''
			gl_longest_orf.write(">" + proper_id + lfr + dupl_status + \
					" STOP:" + `stop_count` + " ANY(X):" + `any_count` + \
					" LONGEST:" + `longest_fr` + " TOTAL_LENGTH:" + `total_size` + \
					" " + full_frame + ":" + `frm_number` + '\n')
			#print longest_orf
			gl_longest_orf.write(longest_orf + '\n')
			'''
			gl_longest_orf += ">" + proper_id + lfr + '\n' + longest_orf + '\n'
		## 5
		if llen5 >= llen1 and llen5 >= llen2 and llen5 >= llen3 and llen5 >= llen4 and llen5 >= llen6:
			if llen5 == llen1 or llen5 == llen2 or llen5 == llen3 or llen5 == llen4 or llen5 == llen6:
				dupl_status = " DUPL_LONG"
			lfr = "_fr5"
			lfrm5 = "fr_5"
			longest_orf = longest_orf5
			stop_count = stop_count5
			any_count = any_count5
			longest_fr = longest_fr5
			total_size = total_size5
			full_frame = full_frame5
			frm_number = frm_number5
			longest_len = llen5
			'''
			gl_longest_orf.write(">" + proper_id + lfr + dupl_status + \
					" STOP:" + `stop_count` + " ANY(X):" + `any_count` + \
					" LONGEST:" + `longest_fr` + " TOTAL_LENGTH:" + `total_size` + \
					" " + full_frame + ":" + `frm_number` + '\n')
			#print longest_orf
			gl_longest_orf.write(longest_orf + '\n')
			'''
			gl_longest_orf += ">" + proper_id + lfr + '\n' + longest_orf + '\n'
		## 6
		if llen6 >= llen1 and llen6 >= llen2 and llen6 >= llen3 and llen6 >= llen4 and llen6 >= llen5:
			if llen6 == llen1 or llen6 == llen2 or llen6 == llen3 or llen6 == llen4 or llen6 == llen5:
				dupl_status = " DUPL_LONG"
			lfr = "_fr6"
			lfrm6 = "fr_6"
			longest_orf = longest_orf6
			stop_count = stop_count6
			any_count = any_count6
			longest_fr = longest_fr6
			total_size = total_size6
			full_frame = full_frame6
			frm_number = frm_number6
			longest_len = llen6
			'''
			gl_longest_orf.write(">" + proper_id + lfr + dupl_status + \
					" STOP:" + `stop_count` + " ANY(X):" + `any_count` + \
					" LONGEST:" + `longest_fr` + " TOTAL_LENGTH:" + `total_size` + \
					" " + full_frame + ":" + `frm_number` + '\n')
			#print longest_orf
			gl_longest_orf.write(longest_orf + '\n')
			'''
			gl_longest_orf += ">" + proper_id + lfr + '\n' + longest_orf + '\n'
		
def InputFileReader(gff_name, bwa_edgename, bwa_pathname, input_read, out_name):
    # Reading FGS predictions on input reads
    gff = open(gff_name, 'r')
    readNameMap = {}
    read_Predictions = set()
    for line in gff.readlines():
        if line[0] != '#':
            line = line.rstrip()
            fields = line.split("\t")
            rname = fields[0];
	    #print rname
            read_Predictions.add(rname)
    # Reading bwa mapping of reads to edges file
    bwaedge = open(bwa_edgename,'r')
    EP_predictions = set()
    for line in bwaedge.readlines():
        if(line[0] != '@' and line[0] != 'S'):
            line = line.rstrip()
            fields = line.split("\t")
            edge_read = fields[0];
            flag = fields[2]
            if flag != '*':
		EP_predictions.add(edge_read)
    # Reading bwa mapping of reads to paths file
    bwapath = open(bwa_pathname,'r')
    for line in bwapath.readlines():
        if(line[0] != '@' and line[0] != 'S'):
            line = line.rstrip()
            fields = line.split("\t")
            path_read = fields[0];
            flag = fields[2]
            if flag != '*':
                  EP_predictions.add(path_read)
    if len(read_Predictions) != 0 and len(EP_predictions) != 0:
        getSequenceInfo(read_Predictions, EP_predictions, fasta_read, out_name)
    gff.close()
    bwaedge.close()
    bwapath.close()


def getSequenceInfo(read_Predictions, EP_predictions, input_read, outd):
    global seqs_to_be_translated
    fin = open(input_read, 'r')
    outfiletag = input_read[0:input_read.rfind('.')]
    predictedOutfile = open (outd+'/reads.filter.fasta', 'w')
    sequence = ""
    seqs_to_be_translated = ""
    allPredictedReadSet = read_Predictions.union(EP_predictions)
    line = fin.readline()
    ReadInfo = {}
    while line:
        if line[0] == '>':
            sequence = ""
            read = line[1:]
	    read = read.rstrip()
            if(read in allPredictedReadSet):
                seq = fin.readline()
		seq = seq.rstrip()
		ReadInfo[read] = seq
                out_header = ">"
                out_header += read
                #sequence = out_header + "\n" + seq
		predictedOutfile.write(out_header)
                predictedOutfile.write('\n')
                predictedOutfile.write(seq)
		predictedOutfile.write('\n')
            else:
                seq = fin.readline()
		seq = seq.rstrip()
		ReadInfo[read] = seq
                out_header = ">"
                out_header += read
		sequence = out_header + '\n' + seq + '\n'
        
	    seqs_to_be_translated += sequence
        line = fin.readline()

    if seqs_to_be_translated != "" :
    	sixFrameTranslate(seqs_to_be_translated, outfiletag)
    else:
	gl_longest_orf = ""
 
    
    if gl_longest_orf != "" :
	lines = gl_longest_orf.split('\n')
	for line in lines :
		fasta_match = line[0:1]
		if fasta_match == ">" :
			read_id = line[1:]
			read_id = read_id[0:read_id.find("_")]
			if read_id not in allPredictedReadSet :
				predictedOutfile.write(">"+read_id+'\n'+ReadInfo[read_id]+'\n')
    predictedOutfile.close()
    


def sixFrameTranslate(seqs_to_be_translated, outfiletag):
    	global  gl_longest_orf
	global seq_type, gen_code, trans_fr
        seq_type =  "DNA"
        gen_code = 1
        trans_fr = 6
    	seqs_min_len = 24
	agct_list   = ["A", "G", "C", "T"]
    	binary_list = ["00", "01", "10", "11"]

	if seq_type == "DNA":
		if gen_code == 1:
			CodonTranslation()
            	if trans_fr == 6:
			gl_longest_orf = ""
                	#gl_longest_orf = open(outfiletag + '.longest_frame', "wb")
    	fasta_id_array = []
	line_counter = 0
	have_seqs = ""
	proper_id = ""
	my_seqs = []
	tot_len = 0
	a_tot = 0
	b_tot = 0
	c_tot = 0
	d_tot = 0
	e_tot = 0
	f_tot = 0
	g_tot = 0
	h_tot = 0
	i_tot = 0
	j_tot = 0
	k_tot = 0
	l_tot = 0
	m_tot = 0
	n_tot = 0
	o_tot = 0
	p_tot = 0
	q_tot = 0
	r_tot = 0
	s_tot = 0
	t_tot = 0
	u_tot = 0
	v_tot = 0
	w_tot = 0
	x_tot = 0
	y_tot = 0
	z_tot = 0
	stop_tot = 0
	seqs_to_be_translated = seqs_to_be_translated.rstrip()
    	lines = seqs_to_be_translated.split('\n')
    	if lines != "" :
		for t in lines:
			t = t.rstrip()
			if t == '':
				###  SUB_SEQ FUNCTION  ###
				have_seqs = "".join(my_seqs)
				seqs_len = len(have_seqs)
				###   STRING PROCESSING   ###
				abc_up = ""
				if seqs_len != 0:
					have_seqs = "".join(my_seqs)
					seqs_len = len(have_seqs)
					tot_len = tot_len + seqs_len
					if seqs_len == 0:
						seqs_len = 1
					###   STRING PROCESSING   ###
					if seq_type == "DNA":
						have_seqs = NSeq_Processor(have_seqs)
	                		abc_up = have_seqs.upper()
	                		a_count = abc_up.count("A")
					a_tot = a_tot + a_count
					###
					b_count = abc_up.count("B")
					b_tot = b_tot + b_count
					###
					c_count = abc_up.count("C")
					c_tot = c_tot + c_count
					###
					d_count = abc_up.count("D")
					d_tot = d_tot + d_count
					###
					e_count = abc_up.count("E")
					e_tot = e_tot + e_count
					###
					f_count = abc_up.count("F")
					f_tot = f_tot + f_count
					###
					g_count = abc_up.count("G")
					g_tot = g_tot + g_count
					###
					h_count = abc_up.count("H")
					h_tot = h_tot + h_count
					###
					i_count = abc_up.count("I")
					i_tot = i_tot + i_count
					###
					j_count = abc_up.count("J")
					j_tot = j_tot + j_count
					###
					k_count = abc_up.count("K")
					k_tot = k_tot + k_count
					###
					l_count = abc_up.count("L")
					l_tot = l_tot + l_count
					###
					m_count = abc_up.count("M")
					m_tot = m_tot + m_count
					###
					n_count = abc_up.count("N")
					n_tot = n_tot + n_count
					###
					o_count = abc_up.count("O")
					o_tot = o_tot + o_count
					###
					p_count = abc_up.count("P")
					p_tot = p_tot + p_count
					###
					q_count = abc_up.count("Q")
					q_tot = q_tot + q_count
					###
					r_count = abc_up.count("R")
					r_tot = r_tot + r_count
					###
					s_count = abc_up.count("S")
					s_tot = s_tot + s_count
					###
					t_count = abc_up.count("T")
					t_tot = t_tot + t_count
					###
					u_count = abc_up.count("U")
					u_tot = u_tot + u_count
					###
					v_count = abc_up.count("V")
					v_tot = v_tot + v_count
					###
					w_count = abc_up.count("W")
					w_tot = w_tot + w_count
					###
					x_count = abc_up.count("X")
					x_tot = x_tot + x_count
					###
					y_count = abc_up.count("Y")
					y_tot = y_tot + y_count
					###
					z_count = abc_up.count("Z")
					z_tot = z_tot + z_count
					###
					stop_count = abc_up.count("*")
					stop_tot = stop_tot + stop_count
					###
					a_fract = round((a_count*100.00/seqs_len),2)
					b_fract = round((b_count*100.00/seqs_len),2)
					c_fract = round((c_count*100.00/seqs_len),2)
					d_fract = round((d_count*100.00/seqs_len),2)
					e_fract = round((e_count*100.00/seqs_len),2)
					f_fract = round((f_count*100.00/seqs_len),2)
					g_fract = round((g_count*100.00/seqs_len),2)
					h_fract = round((h_count*100.00/seqs_len),2)
					i_fract = round((i_count*100.00/seqs_len),2)
					j_fract = round((j_count*100.00/seqs_len),2)
					k_fract = round((k_count*100.00/seqs_len),2)
					l_fract = round((l_count*100.00/seqs_len),2)
					m_fract = round((m_count*100.00/seqs_len),2)
					n_fract = round((n_count*100.00/seqs_len),2)
					o_fract = round((o_count*100.00/seqs_len),2)
					p_fract = round((p_count*100.00/seqs_len),2)
					q_fract = round((q_count*100.00/seqs_len),2)
					r_fract = round((r_count*100.00/seqs_len),2)
					s_fract = round((s_count*100.00/seqs_len),2)
					t_fract = round((t_count*100.00/seqs_len),2)
					u_fract = round((u_count*100.00/seqs_len),2)
					v_fract = round((v_count*100.00/seqs_len),2)
					w_fract = round((w_count*100.00/seqs_len),2)
					x_fract = round((x_count*100.00/seqs_len),2)
					y_fract = round((y_count*100.00/seqs_len),2)
					z_fract = round((z_count*100.00/seqs_len),2)
					stop_fract = round((stop_count*100.00/seqs_len),2)
					###
					at_fract = a_fract + t_fract
					gc_fract = g_fract + c_fract
					atgc_fract = a_fract + t_fract + g_fract + c_fract
					# atgc_fract = round(atgc_fract)
					atgc_fract = round(atgc_fract,2)
					### STRING ###
					a_fract = str(a_fract)
					b_fract = str(b_fract)
					c_fract = str(c_fract)
					d_fract = str(d_fract)
					e_fract = str(e_fract)
					f_fract = str(f_fract)
					g_fract = str(g_fract)
					h_fract = str(h_fract)
					i_fract = str(i_fract)
					j_fract = str(j_fract)
					k_fract = str(k_fract)
					l_fract = str(l_fract)
					m_fract = str(m_fract)
					n_fract = str(n_fract)
					o_fract = str(o_fract)
					p_fract = str(p_fract)
					q_fract = str(q_fract)
					r_fract = str(r_fract)
					s_fract = str(s_fract)
					t_fract = str(t_fract)
					u_fract = str(u_fract)
					v_fract = str(v_fract)
					w_fract = str(w_fract)
					x_fract = str(x_fract)
					y_fract = str(y_fract)
					z_fract = str(z_fract)
					stop_fract = str(stop_fract)

					at_fract = str(at_fract)
					gc_fract = str(gc_fract)
					atgc_fract = str(atgc_fract)
	        		        if trans_fr == 6 and seqs_len >= seqs_min_len:
						Seqs_Translator(proper_id, have_seqs, trans_fr, gen_code)
				if have_seqs != "" and seqs_len >= seqs_min_len:
					break

			if '\n' in t:
				t = t[:-1]
			if '\r' in t:
				t = t[:-1]
			fasta_match = t[0:1]
			if fasta_match == ">":
				#print t
				gi_test = t[0:4]
				if gi_test == ">gi|":
					descr_line = t
					descr_line = re.sub('\t', " ", descr_line)
					descr_line = re.sub("^>gi\|", "", descr_line)
					descr_line = re.sub("\|", '\t', descr_line, 1)
	   				line_counter += 1
				else:
					descr_line = t
					descr_line = re.sub('\t', " ", descr_line)
					descr_line = re.sub("^>", "", descr_line)
					descr_line = re.sub("\|", " ", descr_line, 1)
					descr_line = re.sub(" ", '\t', descr_line, 1)
					line_counter = line_counter + 1
				good_head = string.split(descr_line, '\t')[0]
				try:
					long_tail = string.split(descr_line, '\t')[1]
				except:
					long_tail = ""
				dupl_status = "GOOD"
				if good_head in fasta_id_array:
					dupl_status = "BAD"
					running_text = "\n Check input for duplications  \n ID: " + good_head + "\n"
					print running_text
					break
				fasta_id_array.append(good_head)
				if line_counter != 1:
					abc_up = ""
					have_seqs = "".join(my_seqs)
					seqs_len = len(have_seqs)
					tot_len = tot_len + seqs_len
					if seqs_len == 0:
						seqs_len = 1
					if seq_type == "DNA":
						have_seqs = NSeq_Processor(have_seqs)
					abc_up = have_seqs.upper()
					a_count = abc_up.count("A")
					a_tot = a_tot + a_count

					b_count = abc_up.count("B")
					b_tot = b_tot + b_count

					c_count = abc_up.count("C")
					c_tot = c_tot + c_count

					d_count = abc_up.count("D")
					d_tot = d_tot + d_count

					e_count = abc_up.count("E")
					e_tot = e_tot + e_count

					f_count = abc_up.count("F")
					f_tot = f_tot + f_count

					g_count = abc_up.count("G")
					g_tot = g_tot + g_count

					h_count = abc_up.count("H")
					h_tot = h_tot + h_count

					i_count = abc_up.count("I")
					i_tot = i_tot + i_count

					j_count = abc_up.count("J")
					j_tot = j_tot + j_count

					k_count = abc_up.count("K")
					k_tot = k_tot + k_count

					l_count = abc_up.count("L")
					l_tot = l_tot + l_count

					m_count = abc_up.count("M")
					m_tot = m_tot + m_count

					n_count = abc_up.count("N")
					n_tot = n_tot + n_count

					o_count = abc_up.count("O")
					o_tot = o_tot + o_count

					p_count = abc_up.count("P")
					p_tot = p_tot + p_count

					q_count = abc_up.count("Q")
					q_tot = q_tot + q_count

					r_count = abc_up.count("R")
					r_tot = r_tot + r_count

					s_count = abc_up.count("S")
					s_tot = s_tot + s_count

					t_count = abc_up.count("T")
					t_tot = t_tot + t_count

					u_count = abc_up.count("U")
					u_tot = u_tot + u_count

					v_count = abc_up.count("V")
					v_tot = v_tot + v_count

					w_count = abc_up.count("W")
					w_tot = w_tot + w_count

					x_count = abc_up.count("X")
					x_tot = x_tot + x_count

					y_count = abc_up.count("Y")
					y_tot = y_tot + y_count

					z_count = abc_up.count("Z")
					z_tot = z_tot + z_count

					stop_count = abc_up.count("*")
					stop_tot = stop_tot + stop_count

					a_fract = round((a_count*100.00/seqs_len),2)
					b_fract = round((b_count*100.00/seqs_len),2)
					c_fract = round((c_count*100.00/seqs_len),2)
					d_fract = round((d_count*100.00/seqs_len),2)
					e_fract = round((e_count*100.00/seqs_len),2)
					f_fract = round((f_count*100.00/seqs_len),2)
					g_fract = round((g_count*100.00/seqs_len),2)
					h_fract = round((h_count*100.00/seqs_len),2)
					i_fract = round((i_count*100.00/seqs_len),2)
					j_fract = round((j_count*100.00/seqs_len),2)
					k_fract = round((k_count*100.00/seqs_len),2)
					l_fract = round((l_count*100.00/seqs_len),2)
					m_fract = round((m_count*100.00/seqs_len),2)
					n_fract = round((n_count*100.00/seqs_len),2)
					o_fract = round((o_count*100.00/seqs_len),2)
					p_fract = round((p_count*100.00/seqs_len),2)
					q_fract = round((q_count*100.00/seqs_len),2)
					r_fract = round((r_count*100.00/seqs_len),2)
					s_fract = round((s_count*100.00/seqs_len),2)
					t_fract = round((t_count*100.00/seqs_len),2)
					u_fract = round((u_count*100.00/seqs_len),2)
					v_fract = round((v_count*100.00/seqs_len),2)
					w_fract = round((w_count*100.00/seqs_len),2)
					x_fract = round((x_count*100.00/seqs_len),2)
					y_fract = round((y_count*100.00/seqs_len),2)
					z_fract = round((z_count*100.00/seqs_len),2)
					stop_fract = round((stop_count*100.00/seqs_len),2)

					at_fract = a_fract + t_fract
					gc_fract = g_fract + c_fract
					atgc_fract = a_fract + t_fract + g_fract + c_fract

					atgc_fract = round(atgc_fract,2)

					a_fract = str(a_fract)
					b_fract = str(b_fract)
					c_fract = str(c_fract)
					d_fract = str(d_fract)
					e_fract = str(e_fract)
					f_fract = str(f_fract)
					g_fract = str(g_fract)
					h_fract = str(h_fract)
					i_fract = str(i_fract)
					j_fract = str(j_fract)
					k_fract = str(k_fract)
					l_fract = str(l_fract)
					m_fract = str(m_fract)
					n_fract = str(n_fract)
					o_fract = str(o_fract)
					p_fract = str(p_fract)
					q_fract = str(q_fract)
					r_fract = str(r_fract)
					s_fract = str(s_fract)
					t_fract = str(t_fract)
					u_fract = str(u_fract)
					v_fract = str(v_fract)
					w_fract = str(w_fract)
					x_fract = str(x_fract)
					y_fract = str(y_fract)
					z_fract = str(z_fract)
					stop_fract = str(stop_fract)

					at_fract = str(at_fract)
					gc_fract = str(gc_fract)
					atgc_fract = str(atgc_fract)
					if trans_fr == 6 and seqs_len >= seqs_min_len:
						Seqs_Translator(proper_id, have_seqs, trans_fr, gen_code)
				have_seqs = ""
				my_seqs = []
			if fasta_match != ">" and fasta_match != "" and dupl_status == "GOOD":
				proper_id = good_head
				good_name = long_tail
				my_seqs.append(t)

	a_fract_tot = round((a_tot*100.00/tot_len),2)
	b_fract_tot = round((b_tot*100.00/tot_len),2)
	c_fract_tot = round((c_tot*100.00/tot_len),2)
	d_fract_tot = round((d_tot*100.00/tot_len),2)
	e_fract_tot = round((e_tot*100.00/tot_len),2)
	f_fract_tot = round((f_tot*100.00/tot_len),2)
	g_fract_tot = round((g_tot*100.00/tot_len),2)
	h_fract_tot = round((h_tot*100.00/tot_len),2)
	i_fract_tot = round((i_tot*100.00/tot_len),2)
	j_fract_tot = round((j_tot*100.00/tot_len),2)
	k_fract_tot = round((k_tot*100.00/tot_len),2)
	l_fract_tot = round((l_tot*100.00/tot_len),2)
	m_fract_tot = round((m_tot*100.00/tot_len),2)
	n_fract_tot = round((n_tot*100.00/tot_len),2)
	o_fract_tot = round((o_tot*100.00/tot_len),2)
	p_fract_tot = round((p_tot*100.00/tot_len),2)
	q_fract_tot = round((q_tot*100.00/tot_len),2)
	r_fract_tot = round((r_tot*100.00/tot_len),2)
	s_fract_tot = round((s_tot*100.00/tot_len),2)
	t_fract_tot = round((t_tot*100.00/tot_len),2)
	u_fract_tot = round((u_tot*100.00/tot_len),2)
	v_fract_tot = round((v_tot*100.00/tot_len),2)
	w_fract_tot = round((w_tot*100.00/tot_len),2)
	x_fract_tot = round((x_tot*100.00/tot_len),2)
	y_fract_tot = round((y_tot*100.00/tot_len),2)
	z_fract_tot = round((z_tot*100.00/tot_len),2)
	stop_fract_tot = round((stop_tot*100.00/tot_len),2)

    	###
	at_fract_tot = a_fract_tot + t_fract_tot
	gc_fract_tot = g_fract_tot + c_fract_tot
	atgc_fract_tot = a_fract_tot + t_fract_tot + g_fract_tot + c_fract_tot
	# atgc_fract_tot = round(atgc_fract_tot)
	atgc_fract_tot = round(atgc_fract_tot,2)
	### STRING ###
	a_fract_tot = str(a_fract_tot)
	b_fract_tot = str(b_fract_tot)
	c_fract_tot = str(c_fract_tot)
	d_fract_tot = str(d_fract_tot)
	e_fract_tot = str(e_fract_tot)
	f_fract_tot = str(f_fract_tot)
	g_fract_tot = str(g_fract_tot)
	h_fract_tot = str(h_fract_tot)
	i_fract_tot = str(i_fract_tot)
	j_fract_tot = str(j_fract_tot)
	k_fract_tot = str(k_fract_tot)
	l_fract_tot = str(l_fract_tot)
	m_fract_tot = str(m_fract_tot)
	n_fract_tot = str(n_fract_tot)
	o_fract_tot = str(o_fract_tot)
	p_fract_tot = str(p_fract_tot)
	q_fract_tot = str(q_fract_tot)
	r_fract_tot = str(r_fract_tot)
	s_fract_tot = str(s_fract_tot)
	t_fract_tot = str(t_fract_tot)
	u_fract_tot = str(u_fract_tot)
	v_fract_tot = str(v_fract_tot)
	w_fract_tot = str(w_fract_tot)
	x_fract_tot = str(x_fract_tot)
	y_fract_tot = str(y_fract_tot)
	z_fract_tot = str(z_fract_tot)
	stop_fract_tot = str(stop_fract_tot)

	at_fract_tot = str(at_fract_tot)
	gc_fract_tot = str(gc_fract_tot)
	atgc_fract_tot = str(atgc_fract_tot)
	###


    	#print "-------------------------------------------------"
	#print "Completed running six frame translation filter"
	#print "-------------------------------------------------"


if __name__ == "__main__":
    if len(sys.argv) == 6:
		gff_name  = sys.argv[1]
		bwa_edgename = sys.argv[2]
		bwa_pathname = sys.argv[3]
		fasta_read = sys.argv[4]
        	out_name = sys.argv[5]

		InputFileReader(gff_name, bwa_edgename, bwa_pathname, fasta_read, out_name)

