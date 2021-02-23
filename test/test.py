import numpy as np
import matplotlib.pyplot as plt
import plates

def testPlate():
    plate = plates.Plate('dm-d1,2.CSV')
    print(plate.meta)
    print(plate.df)
    print(plate.row('A'))
    print(plate.col(1))

def testUv384():
    plate = plates.UV384m1('dm-d1,2.CSV')
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



def main():
    testPlate()
    testUv384()
if __name__ == '__main__':
    main()
