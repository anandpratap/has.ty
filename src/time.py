from pylab import *
import subprocess
import re
CALL_CMD = "./has.ty --ndata=%d --mdata=%d --nthreads=%d --ncluster=%d --cpu_ifgt --gpu_ifgt | grep 'time'"


def extract_data(ndata, out):
	lines = out.split("\n")
	times = []
	for line in lines:
		t = re.findall(r'\b\d+\b', line)
		if t != []:
			times.append(float(t[0]))
	t_cpu_cluster = times[0]
	t_cpu_ifgt = times[1]
	t_gpu_cluster = times[2]
	t_gpu_source = times[3]
	t_gpu_eval = times[4]
	print t_gpu_eval, t_gpu_source
	speedup = t_cpu_ifgt/(t_gpu_eval + t_gpu_source + 1e-6)
	f = open("time.log", 'a+')
	f.write("%d %f %f %f %f %f %f\n"%(ndata, t_cpu_cluster, t_cpu_ifgt, t_gpu_cluster, t_gpu_source, t_gpu_eval, speedup))
	f.close()

def run(ndata, mdata, nthreads, ncluster):
	proc = subprocess.Popen([CALL_CMD%(ndata, mdata, nthreads, ncluster)], stdout=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print out
	extract_data(ndata, out)

if __name__ == "__main__":
	e = array(range(3, 8))
	for i in e:
		run(int(10**i), int(10**i), 128, 2)

