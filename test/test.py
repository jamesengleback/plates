import numpy as np
import matplotlib.pyplot as plt
import plates

def testPlate():
    plate = plates.Plate('data/dm-d1,2.CSV')
    print(plate.meta)
    print(plate.df)
    print(plate.row('A'))
    print(plate.col(1))


def testUv384():
    plate = plates.UV384m1('data/dm-d1,2.CSV', smiles = ['c1ccccc1'] * 24)
    concs = np.linspace(0,500,8)
    print(plate.blocks)
    data = plate.block(1)
    print(data)
    norm = plates.transforms.normalize(*data, smooth = True)
    print(norm)
    diff = plates.transforms.difference(norm)
    print(diff)
    response = plates.transforms.response(diff)
    print(response)
    km, vmax, c, rsq = plates.anal.MichaelisMenten(concs, response)
    print(km, vmax, c, rsq)
    '''
    fig, ax = plt.subplots(3,2)
    plates.subplotTraces(norm, ax[0][0], concs = concs)
    plates.subplotTraces(diff, ax[0][1], concs = concs)
    plates.subplotMichaelisMenten(concs, response, ax[1][0], km, vmax, c, rsq)
    plates.subplotText(ax[1][1], {'km':round(km, 3), 'vmax':round(vmax, 3), 'rsq':rsq})
    plates.subplotSMILES(ax[2][0], 'c1ccccc1')

    ax[2][1].axis('off')
    plt.tight_layout()
    plt.show()
    '''
    df = plates.reportPlate(plate, 'report-test')
    print(df)

def test_ascii():
    plate = plates.UV384m2('data/3march-echo-wt-assay-ascii.csv', parser = 'ascii')
    data = plate.block(1)
    print(data)

def test_plate_blocks():
    m1 = plates.UV384m1('data/dm-d1,2.CSV')
    m2_ascii = plates.UV384m2('data/3march-echo-wt-assay-ascii.csv', parser = 'ascii')

    m3 = plates.UV384m3('data/3march-echo-wt-assay-ascii.csv', parser = 'ascii')
    print(m1.block(1))
    print(m2_ascii.block(1))

def test_report():
    m1 = plates.UV384m1('data/dm-d1,2.CSV')
    m2= plates.UV384m2('data/3march-echo-wt-assay-ascii.csv', parser = 'ascii')

    m3 = plates.UV384m3('data/3march-echo-wt-assay-ascii.csv', parser = 'ascii')
    m3.report(controlPlate=m2)


def main():
    testPlate()
    #testUv384()
    #test_ascii()
    #test_plate_blocks()
    #test_report()

if __name__ == '__main__':
    main()
