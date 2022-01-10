import Matrix
import math
import matplotlib.pyplot as plt


def u_begin(t):
    return math.sin(t)


def u_end(t):
    return -math.sin(t)


def u_0(x):
    return 0


def solve_f(x, t):
    return math.sin(t) * math.cos(x)


def f(x, t):
    return math.cos(x) * (math.cos(t) + math.sin(t))


class ParabolicSolver:

    def __init__(self, t_beg, t_end, x_beg, x_end, N, K):
        self.t_beg = t_beg
        self.t_end = t_end
        self.x_beg = x_beg
        self.x_end = x_end
        self.h = (x_end - x_beg) / N
        self.r = (t_end - t_beg) / K
        self.x = []
        self.t = []
        self.N = N
        self.K = K
        self.u = [[0 for column in range(N + 1)] for row in range(K + 1)]
        for j in range(0, N + 1):
            self.u[0][j] = u_0(x_beg + j * self.h)


    def explicit(self):
        sigma = self.r / (self.h ** 2)
        for k in range(1, self.K + 1):
            for j in range(1, self.N):
                self.u[k][j] = sigma * self.u[k - 1][j + 1] + (1 - 2 * sigma) * self.u[k - 1][j] + sigma * self.u[k - 1][j - 1] + self.r * f(self.h * j, self.r * k)
            self.u[k][0] = u_begin(k * self.r)
            self.u[k][self.N] = self.h * u_end(k * self.r) + self.u[k][self.N - 1]



    def implicit(self):
        sigma = self.r / (self.h ** 2)
        for k in range(1, self.K + 1):
            m = [[0 for j in range(0, self.N + 1)] for i in range(0, self.N + 1)]
            d = [0 for j in range(0, self.N + 1)]
            m[0][0] = 1
            m[0][1] = 0
            d[0] = u_begin(self.r * k)
            m[self.N][self.N - 1] = -1 / self.h
            m[self.N][self.N] = 1/ self.h
            d[self.N] = u_end(self.r * k)
            for j in range(1, self.N):
                m[j][j - 1] = sigma
                m[j][j] = -(1 + 2 * sigma)
                m[j][j + 1] = sigma
                d[j] = -self.u[k - 1][j] - self.r * f(j * self.h, k * self.r)
            matrix = Matrix.Matrix(l=m)
            row = matrix.tridiag_solve_SLAU(d)
            self.u[k] = row

    def crank_nicholson_solve(self):
        sigma = self.r / (self.h ** 2)
        for k in range(1, self.K + 1):
            m = [[0 for j in range(0, self.N + 1)] for i in range(0, self.N + 1)]
            d = [0 for j in range(0, self.N + 1)]
            m[0][0] = 1
            m[0][1] = 0
            d[0] = u_begin(self.r * k)
            m[self.N][self.N - 1] = -1 / self.h
            m[self.N][self.N] = 1 / self.h
            d[self.N] = u_end(self.r * k)
            for j in range(1, self.N):
                m[j][j - 1] = 1/2 * sigma
                m[j][j] = -(1 + sigma)
                m[j][j + 1] = 1/2 * sigma
                d[j] = -self.u[k - 1][j] - self.r * f(j * self.h, k * self.r) - 1/2 * sigma * (self.u[k - 1][j + 1] - 2 * self.u[k - 1][j] + self.u[k - 1][j - 1])
            matrix = Matrix.Matrix(l=m)
            row = matrix.tridiag_solve_SLAU(d)
            self.u[k] = row

    def plot(self, x_numb=None, t_numb=None):
        if t_numb is None and x_numb is not None:
            ans = [solve_f(self.h * x_numb, i * self.r) for i in range(self.K + 1)]
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set(ylim=[-1.5, 1.5], xlim=[0, 10])
            ax.set_title('')
            ax.plot([0, 0], [-1000, 1000], color='black')
            ax.plot([-1000, 1000], [0, 0], color='black')
            ax.plot([self.r * i for i in range(self.K + 1)], ans, color='black')
            column = [self.u[i][x_numb] for i in range(self.K + 1)]
            ax.plot([self.r * i for i in range(self.K + 1)], column, color='green')
            plt.grid(True)
            error = 0
            for i in range(self.K + 1):
                error += (abs(ans[i] - column[i])) ** 2
            print(math.sqrt(error))
        if x_numb is None and t_numb is not None:
            ans = [solve_f(self.h * j, self.r * t_numb) for j in range(self.N + 1)]
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set(ylim=[-1.5, 1.5], xlim=[0, math.pi/2])
            ax.set_title('')
            ax.plot([0, 0], [-1000, 1000], color='black')
            ax.plot([-1000, 1000], [0, 0], color='black')
            ax.plot([self.h * i for i in range(self.N + 1)], ans, color='black')
            column = [self.u[t_numb][i] for i in range(self.N + 1)]
            ax.plot([self.h * i for i in range(self.N + 1)], column, color='green')
            plt.grid(True)
            error = 0
            for i in range(self.N + 1):
                error += abs(ans[i] - column[i]) ** 2
            print(math.sqrt(error))

#






a = ParabolicSolver(0, 10, 0, math.pi / 2, 10, 100)
a.implicit()
a.plot(t_numb=4)
plt.show()








