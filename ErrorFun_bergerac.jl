# Writing function seeding_error()
# input arguements: seeding_start, seeding_end, seed_date
# seeding_start and seeding_end are in date format "yyyy-mm-dd".
# computing set of Errors:
# for daily cases      : R2, RMSE, NSE, and KGE12

function seeding_error(seeding_start::Date, seeding_end::Date)
    # Ensure seeding_start is before or equal to seeding_end
    if seeding_start > seeding_end
        throw(ArgumentError("seeding_start should be less than or equal to seeding_end"))
    end

    # Initialize an empty vector to store the rmse
    # daily
    date_rmse  = []
    date_NSE   = []
    date_kge   = []

    # computing R2
    date_R2   = []

    # Loop from seeding_start to seeding_end inclusive
    for date in seeding_start:seeding_end
        # simulate model
        # Define start date as a string
        start_date_str = "18.11.2021"
        seed_date_str = string(Dates.format(Date(date), "dd.mm.yyyy"))
        sol1 = chikungunya(44.87, 0.49, start_date_str, "30.09.2025", seed_date_str)
        
        # Get daily cases
        # crop simulation output to comparable range
        sol1_crop = to(from(sol1, Date(2025, 05, 12)), Date(2025, 10, 12))
        sol1_crop_week = collapse(sol1_crop, week, last)
        # get the weekly case numbers
        weekly_sim = diff(sol1_crop_week["i_inst"] .+ sol1_crop_week["Hr"])
        
        # city data
        timestamps = [Date("2025-06-29", "yyyy-mm-dd"), Date("2025-07-06", "yyyy-mm-dd"), Date("2025-07-13", "yyyy-mm-dd"),
                      Date("2025-07-20", "yyyy-mm-dd"), Date("2025-07-27", "yyyy-mm-dd"), Date("2025-08-03", "yyyy-mm-dd"),
                      Date("2025-08-10", "yyyy-mm-dd")]#, Date("2025-08-17", "yyyy-mm-dd")]
        data = [1, 0, 1, 0, 1, 7, 5]#, 19]
        city_ts = TimeArray(timestamps, data, ["cases"])

        # Extract the timestamps from city_ts
        matching_dates = timestamp(city_ts)
        # Select only those entries in daily_case that match the timestamps in city_ts
        city_sim = weekly_sim[matching_dates]

        # data
        dx = values(city_ts["cases"])
        dy = values(city_sim["i_inst_Hr"])

        # compute NSE
        derror = dx .- dy
        dvar   = dx .- mean(dx)

        NSE    = 1 - (sum(derror.^2)/ sum(dvar.^2))
        
        # compute KGE12
        
        # (Pearson) correlation coefficient
        r = cor(dx, dy)                                # pearson correlation coefficient
        beta =  mean(dy)/mean(dx)                      # Bias  
        gamma = (std(dy)/mean(dy))/(std(dx)/mean(dx))  # variability

        KGE12 = 1 - sqrt((r-1)^2 + (beta-1)^2 + (gamma-1)^2)
        
        # R2
        R2 = r^2

        # RMSE
        RMSE  = sqrt(mean(derror.^2))

        push!(date_rmse, RMSE)
        push!(date_NSE, NSE)
        push!(date_kge, KGE12)
        push!(date_R2, R2)
    end

    # Error data to TimeArray
    error_data = (datetime = seeding_start:seeding_end,
    R2    = date_R2,
    RMSE  = date_rmse,
    NSE   = date_NSE,
    KGE12 = date_kge)
    error_ta = TimeArray(error_data; timestamp = :datetime)

    return error_ta
end