@echo %DATE% %TIME%
@SET COUNT=0
:Loop
	@SET /A COUNT+=1
	@echo %COUNT% %DATE% %TIME%
    @java -Xmx500M -cp .;swarmops.jar Stairway_plot_theta_estimation02 ./testdata/two-epoch%COUNT% 1 5000
	@IF "%COUNT%" == "200" GOTO End
    @GOTO Loop
:End
@echo %DATE% %TIME%