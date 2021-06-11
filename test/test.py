import plates

def testUv384m1(path):
    plate = plates.UV384m1(path, concs = [1,2,3,5,6,7,8])
    x = plate.blocks
    #print(x[1].df[1])

def main():
    testUv384m1('data/dm-d1,2.CSV')

if __name__ == '__main__':
    main()
