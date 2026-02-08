import pandas as pd
import matplotlib.pyplot as plt

df1 = pd.read_csv("graph1.csv")

plt.figure()
plt.plot(df1["t"], df1["F0"], label="F0 (точное)")
plt.plot(df1["t"], df1["N"],  label="N (HyperLogLog)")
plt.xlabel("t (доля обработанного потока)")
plt.ylabel("количество уникальных")
plt.title("График №1: точное F0(t) и оценка N(t)")
plt.legend()
plt.grid(True)
plt.savefig("graph1.png", dpi=200)

df2 = pd.read_csv("graph2.csv")

plt.figure()
plt.plot(df2["t"], df2["meanN"], label="E[N(t)]")
plt.fill_between(df2["t"], df2["lower"], df2["upper"], alpha=0.3, label="E[N(t)] ± σ")
plt.xlabel("t (доля обработанного потока)")
plt.ylabel("количество уникальных")
plt.title("График №2: среднее и разброс оценки HyperLogLog")
plt.legend()
plt.grid(True)
plt.savefig("graph2.png", dpi=200)

plt.show()