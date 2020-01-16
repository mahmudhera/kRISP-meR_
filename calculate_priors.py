import numpy as np
import math
import pandas as pd

threshold = 0.01

class poisson:
    def __init__(self, mean):
        self.mean = mean
    def get_probability(self, k):
        total = -self.mean + k*np.log(self.mean)
        total -= sum([math.log(i) for i in range(2, k+1)])
        return np.exp(total)
    def update_mean(self, mean):
        self.mean = mean
    def get_mean(self):
        return self.mean
    
class geometric:
    def __init__(self, mean):
        self.mean = mean
    def get_probability(self, k):
        p = 1.0/self.mean
        return p*(1-p)**(k-1)
    def update_mean(self, mean):
        self.mean = mean
    def get_mean(self):
        return self.mean
    
def read_histogram(filename):
    histo_data = pd.read_csv(filename, delimiter=' ', names = ['index', 'value'])
    dic = {}
    for t in zip(np.array(histo_data['index']), np.array(histo_data['value'])):
        dic[t[0]]=t[1]
    return dic

def determine_points(histo_data):
    lower = higher = -1
    for i in range(2,max(histo_data.keys())):
        if histo_data[i-1] > histo_data[i] and histo_data[i+1] > histo_data[i]:
            lower = i
        if histo_data[i-1] < histo_data[i] and histo_data[i+1] < histo_data[i]:
            higher = i
            break
    return lower, higher

def determine_inversion_point(histo_data):
    return determine_points(histo_data)[0]

def determine_kmer_coverage(histo_data):
    return determine_points(histo_data)[1]

def expectation_maximization(inv_point, data_points, init_mean, max_count=50):
    data_points = {k:v for k,v in data_points.items() if k <= (max_count+1)*init_mean and k >= inv_point}
    poissons = []
    for x in range(1, max_count+1):
        poissons.append(poisson(float(x*init_mean)))
    probabilities = [{k:0 for k in data_points.keys()} for i in range(max_count)]
    tendencies = [{k:0 for k in data_points.keys()} for i in range(max_count)]
    priors = [1.0/max_count for i in range(max_count)]
    iteration_count = 1
    previous_lambda = init_mean
    previous_error_mean = poissons[0].get_mean()
    while True:
        print('Running iteration ' + str(iteration_count))
        print('Calculating new probabilities.')
        for j in data_points.keys():
            for i in range(max_count):
                probabilities[i][j] = poissons[i].get_probability(j)
        print('Calculating tendencies to belong to the poissons.')
        for i in range(max_count):
            for j in data_points.keys():
                total = 0.0
                for ii in range(max_count):
                    total = total + probabilities[ii][j]*priors[ii]
                tendencies[i][j] = 1.0*probabilities[i][j]*priors[i]/total
        print('Calculating the new means.')
        for i in range(max_count):
            total1 = 0.0
            total2 = 0.0
            for j in data_points.keys():
                total1 = total1 + tendencies[i][j]*j*data_points[j]
                total2 = total2 + tendencies[i][j]*data_points[j]
            poissons[i].update_mean(1.0*total1/total2)
        iteration_count = iteration_count + 1
        estimated_averaged_lambda = poissons[0].get_mean()
        for i in range(max_count):
            poissons[i].update_mean((i+1)*estimated_averaged_lambda)
        val1 = 0.0
        val2 = 0.0
        for i in range(0,max_count):
            for j in data_points.keys():
                val1 += tendencies[i][j]*data_points[j]
                val2 += data_points[j]
            priors[i] = val1/val2
        print('Finished iteration ' + str(iteration_count))
        print('Coverage at this iteration: ' + str(estimated_averaged_lambda))
        if abs(previous_lambda-estimated_averaged_lambda)<threshold:
            break
        else:
            previous_lambda = estimated_averaged_lambda
            previous_error_mean = poissons[0].get_mean()
    return poissons[0].get_mean()

def calculate_priors_with_em(inv_point, data, initial_mean, max_count=50):
    coverage = expectation_maximization(inv_point, data, initial_mean, max_count)
    data_points = {k:v for k,v in data.items() if k <= max_count*initial_mean}
    poissons = [geometric(1.01)]
    for x in range(1,max_count):
        poissons.append(poisson(float(x*coverage)))
    priors = []
    for i in range(max_count):
        sum1 = 0
        sum2 = 0
        for j in data_points.keys():
            sum2 += data_points[j]
            if i == 0 and j < inv_point:
                sum1 += data_points[j]*poissons[i].get_probability(j)
            if i != 0 and j >= inv_point:
                sum1 += data_points[j]*poissons[i].get_probability(j)
        priors.append(sum1/sum2)
    return poissons, [i/sum(priors) for i in priors]

def calculate_priors_without_em(inv_point, data, initial_mean, max_count=50):
    print ('finding priors with preset initial coverage value')
    data_points = {k:v for k,v in data.items() if k <= max_count*initial_mean}
    poissons = [geometric(1.01)]
    for x in range(1,max_count):
        poissons.append(poisson(float(x*initial_mean)))
    priors = []
    for i in range(max_count):
        sum1 = 0
        sum2 = 0
        for j in data_points.keys():
            sum2 += data_points[j]
            if i == 0 and j < inv_point:
                sum1 += data_points[j]*poissons[i].get_probability(j)
            if i != 0 and j >= inv_point:
                sum1 += data_points[j]*poissons[i].get_probability(j)
        priors.append(sum1/sum2)
    return poissons, [i/sum(priors) for i in priors]

def determine_posteriors(max_k, distributions, priors):
    prob_k = [0.0 for k in range(max_k)]
    for k in range(1,max_k):
        for i in range(len(distributions)):
            prob_k[k] += distributions[i].get_probability(k)*priors[i]
    return prob_k

def determine_priors_posteriors(histo_data, perform_em, max_priors=50):
    """
    determines the priors, posteriors and the read-coverage
    :param histo_data: dictionary containing histogram data
    :param perform_em: whether to perform EM or not
    :param max_priors: number of poissons
    :return: priors in an array, posteriors in an array, read-coverage, the inversion point
    """
    #data = read_histogram(histo_filename)
    data = histo_data
    inv_point = determine_inversion_point(data)
    init_coverage = determine_kmer_coverage(data)
    if perform_em:
        distributions, priors = calculate_priors_with_em(inv_point, data, 28, max_priors)
    else:
        distributions, priors = calculate_priors_without_em(inv_point, data, init_coverage, max_priors)
    sum_priors = sum(priors)
    lst = []
    for p in priors:
        lst.append(p/sum_priors)
    priors = lst
    posteriors = determine_posteriors(max_priors*init_coverage, distributions, priors)
    sum_posteriors = sum(posteriors)
    lst = []
    for p in posteriors:
        lst.append(p/sum_posteriors)
    posteriors = lst
    print(posteriors)
    return priors, posteriors, init_coverage, inv_point

if __name__=='__main__':
    priors, posteriors, distributions = determine_priors_posteriors('histo', False, 20)
    for k in range(1,100):
        sum_ = 0.0
        for count in range(1, 20):
            prob = distributions[count].get_probability(k)
            prior = priors[count]
            posterior = posteriors[k]
            sum_ += prob*prior*count/posterior
        print(k, sum_)
