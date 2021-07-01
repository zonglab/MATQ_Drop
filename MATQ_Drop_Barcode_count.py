#!/usr/local/bin/python3.5

import sys

Recovered, barcodeRef = sys.argv[1], sys.argv[2]
HandleRef='TGGTAGGTGGTAGAGA'

barcodeList1 = []
barcodeList2 = []

for line in open(barcodeRef, 'r'):
	line = line.strip().split()
	barcodeList1.append(line[0])
	barcodeList2.append(line[1])

Barcode_count={}

def ten_mismatch(BC,BC_ref):
	temp_count = 0
	for i in range (0,len(BC)):
		if BC[i] == BC_ref[i]:
			temp_count += 1
	if temp_count >= len(BC)-10:
		return True
	else:
		return False

def two_mismatch(BC,BC_ref):
	temp_count = 0
	for i in range (0,len(BC)):
		if BC[i] == BC_ref[i]:
			temp_count += 1
	if temp_count >= len(BC)-2:
		return True
	else:
		return False


def twomismatch_in_list(BC,BC_list):
	Flag = False
	BC_assigned = ''
	for BC_ref in BC_list:
		if two_mismatch(BC, BC_ref):
			Flag = True
			BC_assigned = BC_ref
			break
	return Flag,BC_assigned


Barcode={}
true_to_false_map={}

for line in open(Recovered, 'r'):
	lineS = line.strip().split(' ')
	if len(lineS[0]) == 36:
		bc1, Handle, bc2 = lineS[0][0:10], lineS[0][10:26], lineS[0][26:36]
		if ten_mismatch(Handle,HandleRef):
			Flag1, BC_assigned1 = twomismatch_in_list(bc1, barcodeList1)
			Flag2, BC_assigned2 = twomismatch_in_list(bc2, barcodeList2)
			if Flag1 and Flag2:
				BC_assigned = BC_assigned1 + BC_assigned2
				if BC_assigned not in Barcode:
					Barcode[BC_assigned] = 0
					true_to_false_map[BC_assigned] = []
				Barcode[BC_assigned] += 1
				if (bc1+bc2 != BC_assigned) and ((bc1+bc2) not in true_to_false_map[BC_assigned]):
					true_to_false_map[BC_assigned].append(bc1+bc2)

for key in Barcode:
		corrected_barcodes = ",".join(sorted(true_to_false_map[key]))
		print('%s\t%d\t%s' % (key, Barcode[key], corrected_barcodes))
