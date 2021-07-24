import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as scp

def error(sigma, db, da, a, b, y):
    dy = sigma/2 + 0.001/2
    return ((1/a)*(dy) + (1/a)*(db) + ((y-b)/(a*a))*(da)) 

def function(x, a, b, c):
    return (c + a*np.exp(-x/b))

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

def nearest(a, a0):
    idx = np.abs(a - a0).argmin()
    return idx

def get_halfway_time(time, values, halfway):
    focal = nearest(values, halfway)
    indexes = np.array([focal-3,focal-2,focal-1,focal,focal+1,focal+2,focal+3])
    indexes = indexes[values.size>indexes]
    indexes = indexes[indexes>0]
    index_values = values[indexes]
    index_times = time[indexes]
    ebe = np.polyfit(index_times, index_values, 1, cov=True)
    coof = ebe[0] 
    err = ebe[1]
    time = (halfway-coof[1])/coof[0]
    return time, np.sqrt(np.diag(np.array(err))), coof

def get_hit_time(time, values):
    max_index = np.argmax(values)
    baseline, sigma = get_baseline(time, values)
    tru_val= values[:max_index]
    tru_time= time[:max_index]
    halfway = (values[max_index]+baseline)/2
    time, err, coof = get_halfway_time(tru_time,tru_val, halfway)
    print
    dt = error(sigma, err[1], err[0], coof[0], coof[1], halfway)
    return time, halfway, baseline, dt

def get_baseline(time, values):
    bruh = int(time[0]/(time[0]-time[1]))
    temp = min(abs(time[bruh-1]),abs(time[bruh]),abs(time[bruh+1]))
    if (temp == abs(time[bruh-1])):
        bruh = bruh - 1
    elif (temp == abs(time[bruh + 1])):
        bruh = bruh + 1
    bruh = int(bruh*0.7)
    aaayo = values[:bruh]
    aaayo = reject_outliers(aaayo, 4)
    return np.average(aaayo), np.std(aaayo)

def analze_electrons(time, values):
    max_index = np.argmax(values)
    baseline = values[:int(max_index-60)]
    baseline = values[int(max_index-500):]
    baseline = reject_outliers(baseline)
    avg, sigma = magic_average(baseline, (np.ones(len(baseline))*(0.001)))
    baseline = avg
    left = values[:max_index]
    left_time = time[:max_index]
    halfway = ((values[max_index]+baseline)/2)
    start_time, err , coof= get_halfway_time(left_time, left , halfway)
    right = values[max_index:]
    right_time = time[max_index:]
    stop_time, dump1, dump2 = get_halfway_time(right_time, right , halfway)
    dt = error(sigma, err[1], err[0], coof[0], coof[1], halfway)
    return start_time, stop_time, halfway, baseline, dt
     
def draw_shit(time, counter1, counter2):

    f, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(time, counter2, 'g')
    ax2.plot(time, 300*counter1, 'b')
    start, stop, halfway, baseline_el = analze_electrons(time, counter2)
    hit, halfway_hit, baseline_alpha = get_hit_time(time, counter1) 
    ax1.axvline(x=start)
    ax1.axvline(x=stop)
    ax1.axhline(y=halfway)
    ax1.axhline(y=baseline_el)
    ax2.axhline(y=300*baseline_alpha)
    ax2.axhline(y=300*halfway_hit)
    ax2.axvline(x=hit)
    plt.show()
    plt.clf()

def magic_average(values, errors):
    w = 1/np.multiply(errors, errors)
    up = np.sum(np.multiply(values,w))
    down = np.sum(w)
    return (up/down), np.power(down,-.5)


potentials = []
Ntimes = []
Nerrors =[]
Nspeeds = []
Nspeederr = []
bruh1 = []
bruh2 = []
names = np.array(range(6,31))
names = names*100
name1 = "WN/23_{ayy}/23_{ayy}_{whith}.txt"
name2 = "WN/13_{ayy}/13_{ayy}_{whith}.txt"
for num in names:
    potentials.append(num/(1306))
    Ntime = np.array([])
    Nt13 = np.array([])
    Nt23 = np.array([])
    Tsigmas1 = np.array([])
    Tsigmas2 = np.array([])
    Nerror = np.array([])
    Nspeed = np.array([])
    current_name1 = name1.format(ayy=str(num), whith = "{whith}")
    current_name2 = name2.format(ayy=str(num), whith = "{whith}")
    for i in range(1,21):
        t1 = 0
        t2 = 0
        c11 = 0
        c12 = 0
        c21 = 0
        c22 = 0
        tail = str(int((i-(i%10))/10))+str(i%10)
        curr_name1 = current_name1.format(whith=tail)
        curr_name2 = current_name2.format(whith=tail)
        data1 = pd.read_csv(curr_name1, delimiter = "\t", header = [0,1]).transpose()
        data2 = pd.read_csv(curr_name2, delimiter = "\t", header = [0,1]).transpose()
        if(np.isfinite(np.max(data1.values[2]))):
            t1 = data1.values[0]
            c21 = data1.values[2]
            if(list(list(data1.index)[1])[1] == "(mV)"):
                c11 = (data1.values[1])/1000
            else:
                c11 = data1.values[1]
            start1, stop1, halfway1, baseline_el1, dt1_1 = analze_electrons(t1, c21)
            Tsigmas1 = np.append(Tsigmas1, (stop1-start1))
            hit1, halfway_hit1, baseline_alpha1, dt2_1 = get_hit_time(t1, c11)
            Nt23 = np.append(Nt23,(start1-hit1))
        if(np.isfinite(np.max(data2.values[2]))):
            t2 = data2.values[0]
            c22 = data2.values[2]
            if(list(list(data2.index)[1])[1] == "(mV)"):
                c12 = (data2.values[1])/1000
            else:
                c12 = data2.values[1]
            start2, stop2, halfway2, baseline_el2, dt1_2 = analze_electrons(t2, c22)
            Tsigmas2 = np.append(Tsigmas2, (stop2-start2))
            hit2, halfway_hit2, baseline_alpha2, dt2_2 = get_hit_time(t2, c12)
            Nt13 = np.append(Nt13,(start2-hit2))
                        
    print(len(Nt13))

    
    
    Nt13 = reject_outliers(Nt13, m=1.2)
    print(len(Nt13))
    print("#")
    print(len(Nt23))
    Nt23 = reject_outliers(Nt23, m=1.2)
    print(len(Nt23))
    avg_time, err = (np.average(Nt13)-np.average(Nt23)), (np.std(Nt13)+np.std(Nt23))
    dv = 4.6*err/(avg_time*avg_time)
    avg_speed, speed_err = 4.6/avg_time, dv
    Ntimes.append(avg_time)
    Nerrors.append(err)
    Nspeederr.append(speed_err)
    Nspeeds.append(avg_speed)
    lambda13 = np.multiply(np.multiply(avg_speed,Tsigmas2),np.multiply(avg_speed,Tsigmas2))*(3/(2*9.42))
    lambda13 = reject_outliers(lambda13, m=1.2)
    bruh1.append(np.average(lambda13))
    lambda23 = np.multiply(np.multiply(avg_speed,Tsigmas1),np.multiply(avg_speed,Tsigmas1))*(3/(2*4.82))
    lambda23 = reject_outliers(lambda23, m=1.2)
    bruh2.append(np.average(lambda23))
    #print(np.average(lambda23))
times = Ntimes
errors = Nerrors
dv = np.array(Nspeederr)
plt.errorbar((np.array(potentials)), np.array(times),yerr = np.array(errors))
plt.ylabel("t12 [us]")
plt.xlabel("E/p [V/cm*hPa]")
plt.title("czas przelotu")
print(bruh1)
print(bruh2)
#print(np.array(potentials).tolist())
#print(np.array(times).tolist())
#print(np.array(errors).tolist())
plt.show()
plt.errorbar((np.array(potentials)), Nspeeds,yerr = dv)
plt.ylabel("v12 [cm/us]")
plt.xlabel("E/p [V/cm*hPa]")
plt.title("zależność prędkości od pola elektrycznego")
#print(Nspeeds)
#print(dv.tolist())
plt.show()