import math
import random
from scipy.stats.distributions import chi2
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

seed = 123456789


def linear_congruential_initial(amount):
    a = 101427
    c = 321
    m = 65536
    x_value = seed

    random_numbers = []

    for i in range(0, amount):
        x_value = (a * x_value + c) % m
        random_numbers.append(x_value / m)

    return random_numbers


def linear_congruential_randu(amount):
    a = 65539
    c = 0
    m = 2147483648
    x_value = seed

    random_numbers = []

    for i in range(0, amount):
        x_value = (a * x_value + c) % m
        random_numbers.append(x_value / m)

    return random_numbers


def pythonn_random_numbers(amount):
    random_numbers = []
    for i in range(0, amount):
        random_numbers.append(random.random())

    return random_numbers

# Compares a cdf of the RNs with a cdf of an ideal uniform distribution
def kolmogorov_smirnov_test(numbers):
    r = []
    if len(numbers) > 100:
        for i in range(0, 100):
            r.append(numbers[i])
        else:
            r = numbers

    r.sort()                            
    n = len(r)          

    d_positive = 0
    d_negative = 0

    for i in range(0, len(r)):
        temp_d = (i + 1) / n - r[i]
        if temp_d > d_positive:
            d_positive = temp_d        

        temp_d = (r[i] - i / n)
        if temp_d > d_negative:
            d_negative = temp_d

    d = max(d_positive, d_negative) # Compute the actual distribution of the RNs
    d_alpha = 1.36 / math.sqrt(n) # Ideal uniform distribution tabular value for the siginififance level 0.05 and sample size n

    return d <= d_alpha


def runs_test(numbers):
    actual_run_lenghts = compute_actual_runs(numbers)
    optimal_run_lenghts = compute_optimal_runs(len(numbers))
    chi_squared = 0 
    chi_squared_af = 0
    p = 1
    alpha = 0.05
    for i in range(1, len(actual_run_lenghts) + 1):
        Oi = actual_run_lenghts.get(i, 0)
        Ei = optimal_run_lenghts.get(i, 0)
        if not (Ei == 0 and Oi == 0):
            runs_test = (math.pow(Ei - Oi, 2)) / Ei
            chi_squared += runs_test # Compute the chi-squared value for the RNs
    
    # Use the inverse function of chi-squared to compute the ideal chi-squared tabular value given d.o.f and a
    chi_squared_af = chi2.ppf(p - alpha, len(actual_run_lenghts) - 1)

    # Compare actual chi-sqaured distribution with an ideal chi-squared distribution
    return not (chi_squared > chi_squared_af)


def compute_actual_runs(numbers):
    actual_run_lenghts = {}
    lenght = 0
    ascending = False
    descending = False

    for i in range(0, len(numbers) - 1):
        if numbers[i] > numbers[i + 1] and ascending:
            actual_run_lenghts[lenght] = actual_run_lenghts.get(lenght, 0) + 1
            lenght = 0
            ascending = False
        if numbers[i] < numbers[i + 1] and descending:
            actual_run_lenghts[lenght] = actual_run_lenghts.get(lenght, 0) + 1
            lenght = 0
            descending = False
        if numbers[i] == numbers[i + 1]:
            lenght += 1
        elif numbers[i] > numbers[i + 1]:
            lenght += 1
            descending = True
        elif numbers[i] < numbers[i + 1]:
            lenght += 1
            ascending = True
    actual_run_lenghts[lenght] = actual_run_lenghts.get(lenght, 0) + 1
    return actual_run_lenghts


def compute_optimal_runs(sequenceLength):
    optimal_run_lenghts = {}
    N = sequenceLength
    for i in range(1, N):
        optimalRunLenght = 2 / math.factorial(i + 3) * (
                    N * (math.pow(i, 2) + 3 * i + 1) - (math.pow(i, 3) + 3 * math.pow(i, 2) - i - 4))
        optimal_run_lenghts[i] = optimalRunLenght
    return optimal_run_lenghts


def create_3D_scatterplot(numbers):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    numbers.pop()  # To get 9999 numbers, since 10000 cannot be represented
    x = numbers[0::3]
    y = numbers[1::3]
    z = numbers[2::3]

    ax.scatter(x, y, z, c='r', marker='o')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


standard = pythonn_random_numbers(10000)
# create_3D_scatterplot(standard)

lcm_initial = linear_congruential_initial(10000)
# create_3D_scatterplot(lcm_initial)

lcm_RANDU = linear_congruential_randu(10000)
# create_3D_scatterplot(lcm_RANDU)

print()
print("TESTING FOR UNIFORMITY AND INDEPENDENCE")
print()
print("Standard RNs uniformity: " + str(kolmogorov_smirnov_test(standard)))
print("Standard RNs independence: " + str(runs_test(standard)))
print()
print("Initial RNs uniformity: " + str(kolmogorov_smirnov_test(lcm_initial)))
print("Initial RNs independence: " + str(runs_test(lcm_initial)))
print()
print("Randu RNs uniformity: " + str(kolmogorov_smirnov_test(lcm_RANDU)))
print("Randu RNs independence: " + str(runs_test(lcm_RANDU)))
