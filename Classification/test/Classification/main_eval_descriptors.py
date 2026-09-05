
with open('build/results.txt', mode='r') as f:
    for line in f:
        text = line.strip()
        if '(train) accuracy' in text:
            print('train accuracy:',text.split(',')[1])
        if '(test) accuracy' in text:
            print('test accuracy:',text.split(',')[1])
