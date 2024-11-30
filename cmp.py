#!/bin/python3

fn1 = 'results.log'
fn2 = 't'

f1 = open(fn1, "r").readlines()
f2 = open(fn2, "r").readlines()

i = 2
nb_diff = 0
for l in f1[0:len(f1) - 1]:
    inst_r = l.split()[0]
    bins_r = int(l.split()[1])

    inst_t = f2[i - 1]
    bins_t = int(f2[i].split(':')[1])

    if bins_r != bins_t:
        print(f'Diff insts:{inst_r} {inst_t} bins:{bins_r} {bins_t}')
        nb_diff += 1
    # else:
    #     print(f'Insts {inst_r} {inst_t} Ok')

    i += 5
    if i > len(f2):
        break;

print(f'Nb diff: {nb_diff} of:{len(f1)}')