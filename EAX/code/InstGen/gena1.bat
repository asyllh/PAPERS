@echo off
for /l %%i in (1,1,1000) do (
	echo Creating PA1 - Artficial Instance %%i of 1000 for 100 items
	instgen.exe -i %%i -s %%i -n 100 -c 1
)
echo FINISHED GENERATING PA1 INSTANCES
