#!/usr/bin/python

KEYS='WIDTH HEIGHT ITERS SAVE_PERIOD AVERAGING_RADIUS ENSEMBLE_SIZE INIT1 INIT2 INIT4 INIT8 INIT16 INIT32 INIT64 INIT128 SOURCE_WIDTH SRC1 SRC2 SRC4 SRC8 SRC16 SRC32 SRC64 SRC128'

import sys, os

KEYS=KEYS.split(' ')
failed=False
for key in KEYS:
	if key not in os.environ:
		failed=True
		sys.stderr.write("ERROR: undefined environment variable:" + str(key) + '\n')
assert not failed



print('// Generated by', sys.argv[0])
print('#define THREADS 8')
print('#define WIDTH', os.environ['WIDTH'])
print('#define HEIGHT', os.environ['HEIGHT'])
print('#define ITERS', os.environ['ITERS'])
print('#define SAVE_PERIOD', os.environ['SAVE_PERIOD'])
print('#define AVERAGING_RADIUS', os.environ['AVERAGING_RADIUS'])
print('#define ENSEMBLE_SIZE', os.environ['ENSEMBLE_SIZE'])
print('#define SOURCE_WIDTH', os.environ['SOURCE_WIDTH'])

print('''
// Field initialization
void init()
{
	fill(0, 0, HEIGHT, WIDTH, ''' + \
	', '.join([os.environ[key] for key in [
		'INIT1', 'INIT2', 'INIT4', 'INIT8', 
		'INIT16', 'INIT32', 'INIT64', 'INIT128'
	]]) + ''');
}

// Sources rules
void sources(int iter)
{
	if (iter==0) {
		fill(0, WIDTH/2-SOURCE_WIDTH/2, HEIGHT, WIDTH/2+SOURCE_WIDTH/2, ''' + \
			', '.join([os.environ[key] for key in [
				'SRC1', 'SRC2', 'SRC4', 'SRC8',
				'SRC16', 'SRC32', 'SRC64', 'SRC128'
			]]) +''');
	}
}
''')
