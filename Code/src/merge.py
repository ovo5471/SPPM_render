# -*- coding: UTF-8 -*-
# 参考代码：https://blog.csdn.net/rocketeerLi/article/details/84929516


import cv2
import sys
import numpy as np
from struct import unpack

# 读取并存储 bmp 文件


class ReadBMPFile:
    def __init__(self, filePath):
        file = open(filePath, "rb")
        # 读取 bmp 文件的文件头    14 字节
        # 0x4d42 对应BM 表示这是Windows支持的位图格式
        self.bfType = unpack("<h", file.read(2))[0]
        self.bfSize = unpack("<i", file.read(4))[0]       # 位图文件大小
        self.bfReserved1 = unpack("<h", file.read(2))[0]  # 保留字段 必须设为 0
        self.bfReserved2 = unpack("<h", file.read(2))[0]  # 保留字段 必须设为 0
        # 偏移量 从文件头到位图数据需偏移多少字节（位图信息头、调色板长度等不是固定的，这时就需要这个参数了）
        self.bfOffBits = unpack("<i", file.read(4))[0]
        # 读取 bmp 文件的位图信息头 40 字节
        self.biSize = unpack("<i", file.read(4))[0]       # 所需要的字节数
        self.biWidth = unpack("<i", file.read(4))[0]      # 图像的宽度 单位 像素
        self.biHeight = unpack("<i", file.read(4))[0]     # 图像的高度 单位 像素
        self.biPlanes = unpack("<h", file.read(2))[0]     # 说明颜色平面数 总设为 1
        self.biBitCount = unpack("<h", file.read(2))[0]   # 说明比特数

        self.biCompression = unpack("<i", file.read(4))[0]  # 图像压缩的数据类型
        self.biSizeImage = unpack("<i", file.read(4))[0]    # 图像大小
        self.biXPelsPerMeter = unpack("<i", file.read(4))[0]  # 水平分辨率
        self.biYPelsPerMeter = unpack("<i", file.read(4))[0]  # 垂直分辨率
        self.biClrUsed = unpack("<i", file.read(4))[0]      # 实际使用的彩色表中的颜色索引数
        self.biClrImportant = unpack("<i", file.read(4))[
            0]  # 对图像显示有重要影响的颜色索引的数目
        self.bmp_data = []

        if self.biBitCount != 24:
            print("输入的图片比特值为 ：" + str(self.biBitCount) + "\t 与程序不匹配")

        for height in range(self.biHeight):
            bmp_data_row = []
            # 四字节填充位检测
            count = 0
            for width in range(self.biWidth):
                bmp_data_row.append([unpack("<B", file.read(1))[0], unpack(
                    "<B", file.read(1))[0], unpack("<B", file.read(1))[0]])
                count = count + 3
            # bmp 四字节对齐原则
            while count % 4 != 0:
                file.read(1)
                count = count + 1
            self.bmp_data.append(bmp_data_row)
        self.bmp_data.reverse()
        file.close()
        # R, G, B 三个通道
        self.R = []
        self.G = []
        self.B = []

        for row in range(self.biHeight):
            R_row = []
            G_row = []
            B_row = []
            for col in range(self.biWidth):
                B_row.append(self.bmp_data[row][col][0])
                G_row.append(self.bmp_data[row][col][1])
                R_row.append(self.bmp_data[row][col][2])
            self.B.append(B_row)
            self.G.append(G_row)
            self.R.append(R_row)


# 命令行传入的文件路径
filePath1 = sys.argv[1]
filePath2 = sys.argv[2]
filePath3 = sys.argv[3]
# 读取 BMP 文件
bmpFile1 = ReadBMPFile(filePath1)
bmpFile2 = ReadBMPFile(filePath2)
bmpFile3 = ReadBMPFile(filePath3)
# R, G, B 三个通道 [0, 255]
R = bmpFile1.R
G = bmpFile2.G
B = bmpFile3.B
# 显示图像
b = np.array(B, dtype=np.uint8)
g = np.array(G, dtype=np.uint8)
r = np.array(R, dtype=np.uint8)
merged = cv2.merge([b, g, r])  # 合并R、G、B分量 默认顺序为 B、G、R
cv2.imshow("Merged", merged)
cv2.imwrite("output/merged.bmp", merged)
