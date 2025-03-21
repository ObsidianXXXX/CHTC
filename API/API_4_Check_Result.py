import os
dis_type = "2D Exponential"
BatchDir = os.path.abspath("..\\" + dis_type + "\\Batch\\")

Batch_Num = 20
succ = []
fail = []
for i in range(Batch_Num):
    abqWorkDir = os.path.join(BatchDir, str(i + 1))
    os.chdir(abqWorkDir)
    stress_bin = os.path.join(abqWorkDir, str(i + 1) + 'Stress.bin')
    strain_bin = os.path.join(abqWorkDir, str(i + 1) + 'Strain.bin')
    if os.path.exists(stress_bin) and os.path.exists(strain_bin):
        succ.append(i + 1)
        print(str(i + 1) + ': success')
    else:
        fail.append(i + 1)
        print(str(i + 1) + ': fail')
os.chdir("..\\")
succ_file = os.path.abspath("..\\" + "succ.txt")
fail_file = os.path.abspath("..\\" + "fail.txt")

with open(succ_file,"w") as f:
    for i in succ:
        f.write(str(i) + "\n")

with open(fail_file,"w") as f:
    for i in fail:
        f.write(str(i) + "\n")
        
print("succ:",succ)
print("fail:",fail)
