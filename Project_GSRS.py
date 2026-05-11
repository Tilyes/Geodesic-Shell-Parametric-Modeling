# -*- coding: utf-8 -*-
import math
import copy
import numpy as np
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
"""简单坐标运算"""
def plus(x,y):
	return (x[0] + y[0], x[1] + y[1], x[2] + y[2])
def cut(x,y):
	return (x[0] - y[0], x[1] - y[1], x[2] - y[2])

def double(x,m):
	return (x[0]*m, x[1]*m, x[2]*m)

def distance(x,y):
	ai = x[0] - y[0]
	bi = x[1] - y[1]
	ci = x[2] - y[2]
	d = math.pow((math.pow(ai,2)+math.pow(bi,2)+math.pow(ci,2)), 0.5)
	return d

# 点映射到球面上的点坐标：等价于 r * x / |x|
def xyz_coordinates(x, radius):
	length = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
	if length == 0:
		return (0.0, 0.0, 0.0)
	scale = float(radius) / length
	return (x[0] * scale, x[1] * scale, x[2] * scale)


"""粗略寻找较好的调整系数"""
#运用列举法以求杆长长度的方差最小时的调整系数
def search_best_adc(f1,f2,s,h):

	Kn_number = 5  # 圆分数
	pi = math.pi
	angle = float(2 * pi) / Kn_number
	radius = float(h) / 2 + (math.pow(s, 2) / 8) / h  # 半径
	z_angle = math.acos((radius - h) / radius)
	coefficient = float(f1) / (f1 + f2)  # 上下部分角度的分配系数
	ad_list = [-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12]
	t = radius*radius
	for m in range(len(ad_list)):
		adjustment_coefficient = ad_list[m]
		Distances = []

		last_coefficient = coefficient + adjustment_coefficient  # 最终分配系数
		z_angle_1 = last_coefficient * z_angle
		"""上半部分坐标计算"""
		coordinates0 = [ (math.sin(z_angle_1)*math.cos(angle*m),math.sin(z_angle_1)*math.sin(angle*m),math.cos(z_angle_1)) for m in range(Kn_number) ]
		coordinates1 = [ ]
		m = 0
		for c in coordinates0:
			a = double(cut(c , (0,0,1)), 1.0/f1)
			coordinates = [ plus((0,0,1),double(a,m)) for m in range(f1+1) ]
			coordinates1.append(coordinates)
			m = m + 1

		coordinates2 = [[(0,0,radius)],]
		coordinates = [ ]
		m=1
		while m <= f1:
			i=0
			while i < Kn_number:
				if i+1 < Kn_number :
					a = double(cut(coordinates1[i][m],coordinates1[i+1][m]),1.0/m)
				else:
					a = double(cut(coordinates1[i][m],coordinates1[0][m]),1.0/m)
				for j in range(m):
					coordinates.append(cut(coordinates1[i][m], double(a,j)))
				i = i + 1
			coordinates2.append(coordinates)
			coordinates = []
			m = m + 1

		for m in range(f1+1):
			if m != 0:
				for q in range(Kn_number*m):
					coordinates2[m][q] = xyz_coordinates(coordinates2[m][q],radius)

		"""上半部分杆长的计算"""
		#同一层杆件的生成
		m = 1
		while m <= f1:

			for q in range(Kn_number*m):
				if q+1 < Kn_number*m:
					Distances.append(distance(coordinates2[m][q],coordinates2[m][q+1]))
				else:
					Distances.append(distance(coordinates2[m][q],coordinates2[m][0]))
			m = m + 1

		#相邻层的杆件生成
		coordinates3 = copy.deepcopy(coordinates2) #列表的深层复制
		m = f1
		while m > 0 :
			q = Kn_number*m -1
			while q >= 0:
				if q % m == 0:
					n = q / m * (m-1)
					Distances.append(distance(coordinates2[m][q],coordinates2[m-1][int(n)]))
					del coordinates3[m][q]
				q = q-1
			for q in range(Kn_number*(m-1)):
				if q+1 < Kn_number*(m-1):
					Distances.append(distance(coordinates3[m][q],coordinates2[m-1][q+1]))
					Distances.append(distance(coordinates3[m][q],coordinates2[m-1][q]))
				else:
					Distances.append(distance(coordinates3[m][q],coordinates2[m-1][0]))
					Distances.append(distance(coordinates3[m][q],coordinates2[m-1][q]))
			m = m - 1

		"""下半部分坐标计算"""
		c1 = []
		for i in range(5*f1):
			if i % f1 ==0:
				c1.append(double(coordinates2[f1][i], 1.0/radius))

		#下界限坐标
		c2 = []
		for i in range(5):
			c2.append(  (math.sin(z_angle)*math.cos(float(f2)/f1*pi/5+i*2*pi/5), math.sin(z_angle)*math.sin(float(f2)/f1*pi/5+i*2*pi/5), math.cos(z_angle)) )
		# 下界限辅助计算坐标
		c21 = []
		for i in range(5):
			c21.append(  (math.sin(z_angle)*math.cos(float(f2)/f1*pi/5+(float(f1-f2))/f1*2*pi/5+2*i*pi/5), math.sin(z_angle)*math.sin(float(f2)/f1*pi/5+(float(f1-f2))/f1*2*pi/5+2*i*pi/5), math.cos(z_angle)) )

		#计算主要部分的节点坐标
		m = 0
		n = 0
		c3 = []
		while n < 5:
			a = double(cut(c2[n],c1[m]), 1.0 / f2)
			c0 = [plus(c1[m], double(a, j)) for j in range(f2 + 1)]
			c3.append(c0)
			if abs(m-n)>2:
				break
			if m > n :
				n = n + 1
				continue
			if m == n:
				m = m + 1
				c2 [n] = c21[n]
			if m == 5:
				m = 0

		#计算下半部分所用网点的坐标
		c4 = []
		c00 = []
		for i in range(f2 + 1):
			for m in range(10):
				if m % 2 ==0 and i != f1:
					a = double(cut(c3[m + 1][i], c3[m][i]), 1.0 / (f1 - i))
					for v in range(int(f1) - i):
						c00.append(plus(c3[m][i], double(a, v)))
				elif m %2 != 0 and i!=0 :
					if m < 9:
						a = double(cut(c3[m + 1][i], c3[m][i]), 1.0 / i)
					else:
						a = double(cut(c3[0][i], c3[m][i]), 1.0 / i)
					for v in range(i):
						c00.append(plus(c3[m][i], double(a, v)))
			c4.append(c00)
			c00 = []

		#点投影到球面上
		for i in range(f2 +1):
			for m in range(5*f1):
				c4[i][m] = xyz_coordinates(c4[i][m],radius)

		#!!!!!!上半部分和下半部分相连处的节点统一，使用上半部分的低层坐标
		for i in range(5*f1):
			c4[0][i] = coordinates2[f1][i]


		"""下半部分的杆长计算"""
		#同一层杆件的生成
		m = 1
		while m <= f2:
			for q in range(5* f1):
				if q + 1 < 5 * f1:
					Distances.append(distance(c4[m][q], c4[m][q + 1]))

				else:
					Distances.append(distance(c4[m][q], c4[m][0]))
			m = m+1
		#不同层杆件的生成
		m = 0
		while m < f2:
			for q in range(5 * f1):
				if q > 0:
					Distances.append(distance(c4[m][q], c4[m + 1][q]))
					Distances.append(distance(c4[m][q], c4[m + 1][q-1]))
				else:

					Distances.append(distance(c4[m][q], c4[m + 1][q]))
					Distances.append(distance(c4[m][q], c4[m + 1][-1]))
			m = m + 1

		"""方差的计算"""
		distance_var = np.var(Distances)
		#不同调整系数方差进行比较，取方差最小的时的调整系数
		if distance_var <= t :
			best_adc = adjustment_coefficient
			t = distance_var

	return best_adc



class Project_dcx():
	def __init__(self, part_name, f1, f2 , s, h):

		"""

		参数设置
		part_name = 'duanchengxian'
		f1 = 5  #上半部分层数
		f2 = 3  #下半部分层数，f2 <= f1/2
		s = 200     #跨度
		h = 120     #矢高
		Kn_number = 5     # 圆分数
		angle = float(2*pi)/Kn_number
		radius = float(h)/2 + (math.pow(s,2)/8)/h    #半径
		z_angle = math.acos((radius-h)/radius)
		coefficient = float(f1)/(f1+f2)      #上下部分角度的分配系数
		adjustment_coefficient        # 修正系数，用于调整整体杆长方差
		last_coefficient = coefficient + adjustment_coefficient #最终分配系数
		z_angle_1 = coefficient*z_angle

		"""
		Kn_number = 5
		angle = float(2*pi)/Kn_number
		radius = float(h)/2 + (math.pow(s,2)/8)/h
		z_angle = math.acos((radius-h)/radius)
		coefficient = float(f1)/(f1+f2)
		adjustment_coefficient = search_best_adc(f1,f2,s,h)
		last_coefficient = coefficient + adjustment_coefficient
		z_angle_1 = last_coefficient * z_angle

		mdb.models['Model-1'].Part(name=part_name,dimensionality=THREE_D,type=DEFORMABLE_BODY)
		p = mdb.models['Model-1'].parts[part_name]

		self.f1 = f1
		self.f2 = f2
		self.p = p
		self.radius = radius
		self.z_angle_1 = z_angle_1
		self.z_angle = z_angle
		self.angle = angle
		self.Kn_number = Kn_number

	def part(self):
		pi = math.pi



		"""上半部分坐标计算"""
		#底层坐标计算
		coordinates0 = [ (math.sin(self.z_angle_1)*math.cos(self.angle*m),math.sin(self.z_angle_1)*math.sin(self.angle*m),cos(self.z_angle_1)) for m in range(self.Kn_number) ]

		coordinates1 = [ ]
		m = 0
		#计算棱边上坐标
		for c in coordinates0:
			a = double(cut(c , (0,0,1)), 1.0/self.f1)

			coordinates = [ plus((0,0,1),double(a,m)) for m in range(self.f1+1) ]
			coordinates1.append(coordinates)
			m = m + 1

		coordinates2 = [[(0,0,self.radius)],]  #手动补上顶点坐标
		coordinates = [ ]
		m=1

		#计算各层坐标
		while m <= self.f1:
			i=0
			while i < self.Kn_number:
				if i+1 < self.Kn_number :
					a = double(cut(coordinates1[i][m],coordinates1[i+1][m]),1.0/m)
				else:
					a = double(cut(coordinates1[i][m],coordinates1[0][m]),1.0/m)
				for j in range(m):
					coordinates.append(cut(coordinates1[i][m], double(a,j)))
				i = i + 1
			coordinates2.append(coordinates)
			coordinates = []
			m = m + 1


		#直角坐标转换为球面坐标
		for m in range(self.f1+1):
			if m != 0:
				for q in range(self.Kn_number*m):
					coordinates2[m][q] = xyz_coordinates(coordinates2[m][q],self.radius)

		"""上半部分单元生成"""
		#同一层点的连接
		m = 1
		while m <= self.f1:

			for q in range(self.Kn_number*m):
				if q+1 < self.Kn_number*m:
					self.p.WirePolyLine(points=(coordinates2[m][q],coordinates2[m][q+1]),meshable=ON)
				else:
					self.p.WirePolyLine(points=(coordinates2[m][q],coordinates2[m][0]),meshable=ON)
			m = m + 1

		#相邻层点的连接
		coordinates3 = copy.deepcopy(coordinates2)

		m = self.f1
		while m > 0 :
			q = self.Kn_number*m -1
			while q >= 0:
				if q % m == 0:
					self.p.WirePolyLine(points=(coordinates2[m][q],coordinates2[m-1][q/(m)*(m-1)]),meshable=ON)
					del coordinates3[m][q]
				q = q-1
			for q in range(self.Kn_number*(m-1)):
				if q+1 < self.Kn_number*(m-1):
					self.p.WirePolyLine(points=(coordinates3[m][q],coordinates2[m-1][q+1]),meshable=ON)
					self.p.WirePolyLine(points=(coordinates3[m][q],coordinates2[m-1][q]),meshable=ON)
				else:
					self.p.WirePolyLine(points=(coordinates3[m][q],coordinates2[m-1][0]),meshable=ON)
					self.p.WirePolyLine(points=(coordinates3[m][q],coordinates2[m-1][q]),meshable=ON)
			m = m - 1


		"""下半部分坐标计算"""
		#确定上界限坐标（原正二十面的二等分点及端点）
		c1 = []
		for i in range(self.Kn_number*self.f1):
			if i % self.f1 ==0:
				c1.append(double(coordinates2[self.f1][i], 1.0/self.radius))

		#下界限坐标
		c2 = []
		for i in range(self.Kn_number):
			c2.append(  (math.sin(self.z_angle)*math.cos(float(self.f2)/self.f1*pi/5+i*2*pi/5), math.sin(self.z_angle)*math.sin(float(self.f2)/self.f1*pi/5+i*2*pi/5), math.cos(self.z_angle)) )
		# 下界限辅助计算坐标
		c21 = []
		for i in range(self.Kn_number):
			c21.append(  (math.sin(self.z_angle)*math.cos(float(self.f2)/self.f1*pi/5+(float(self.f1-self.f2))/self.f1*2*pi/5+2*i*pi/5), math.sin(self.z_angle)*math.sin(float(self.f2)/self.f1*pi/5+(float(self.f1-self.f2))/self.f1*2*pi/5+2*i*pi/5), math.cos(self.z_angle)) )

		#辅助列坐标计算
		m = 0
		n = 0
		c3 = []
		while n < self.Kn_number:
			a = double(cut(c2[n],c1[m]), 1.0 / self.f2)
			c0 = [plus(c1[m], double(a, j)) for j in range(self.f2 + 1)]
			c3.append(c0)
			if abs(m-n)>2:
				break
			if m > n :
				n = n + 1
				continue
			if m == n:
				m = m + 1
				c2 [n] = c21[n]
			if m == self.Kn_number:
				m = 0

		#各层坐标计算
		c4 = []
		c00 = []
		for i in range(self.f2 + 1):
			for m in range(2*self.Kn_number):
				if m % 2 ==0 and i != self.f1:
					a = double(cut(c3[m + 1][i], c3[m][i]), 1.0 / (self.f1 - i))
					for v in range(int(self.f1) - i):
						c00.append(plus(c3[m][i], double(a, v)))
				elif m %2 != 0 and i!=0 :
					if m < 2*self.Kn_number - 1:
						a = double(cut(c3[m + 1][i], c3[m][i]), 1.0 / i)
					else:
						a = double(cut(c3[0][i], c3[m][i]), 1.0 / i)
					for v in range(i):
						c00.append(plus(c3[m][i], double(a, v)))
			c4.append(c00)
			c00 = []

		#弦分法的点投影至球面上
		for i in range(self.f2 +1):
			for m in range(5*self.f1):
				c4[i][m] = xyz_coordinates(c4[i][m],self.radius)

		#连接上半部分和下半部分
		for i in range(self.Kn_number*self.f1):
			c4[0][i] = coordinates2[self.f1][i]

		"""下半部分单元生成"""
		#同一层的节点连接（除去顶层）
		m = 1
		while m <= self.f2:
			for q in range(self.Kn_number* self.f1):
				if q + 1 < self.Kn_number * self.f1:
					self.p.WirePolyLine(points=(c4[m][q], c4[m][q + 1]), meshable=ON)

				else:
					self.p.WirePolyLine(points=(c4[m][q], c4[m][0]), meshable=ON)
			m = m+1

		#相邻层的节点连接
		m = 0
		while m < self.f2:
			for q in range(self.Kn_number * self.f1):
				if q > 0:
					self.p.WirePolyLine(points=(c4[m][q], c4[m + 1][q]), meshable=ON)
					self.p.WirePolyLine(points=(c4[m][q], c4[m + 1][q - 1]), meshable=ON)
				else:
					self.p.WirePolyLine(points=(c4[m][q], c4[m + 1][q]), meshable=ON)
					self.p.WirePolyLine(points=(c4[m][q], c4[m + 1][-1]), meshable=ON)
			m = m + 1

def main(part_name, f1, f2 , s, h):

	pj_dcx = Project_dcx(part_name, f1, f2 , s, h)
	pj_dcx.part()
	print("Completed!")

if __name__ == "__main__":
	main(part_name='duanchengxian_project01', f1=5, f2=2 , s=20, h=9)
