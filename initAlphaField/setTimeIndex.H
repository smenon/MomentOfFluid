    instantList sourceTimes = runTimeSource.times();
    label sourceTimeIndex = runTimeSource.timeIndex();

    if (args.options().found("sourceTime"))
    {
        if ((args.options()["sourceTime"]) == "latestTime")
        {
            sourceTimeIndex = sourceTimes.size() - 1;
        }
        else
        {
            sourceTimeIndex = Time::findClosestTimeIndex
            (
                sourceTimes,
                args.optionRead<scalar>("sourceTime")
            );
        }
    }
    else
    {
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            args.optionRead<scalar>("sourceTime")
        );
    }
