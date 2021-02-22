import plates

def main():
    plate = plates.AssayPlate('dm-d1,2.CSV')
    plate.report('test')

if __name__ == '__main__':
    main()
