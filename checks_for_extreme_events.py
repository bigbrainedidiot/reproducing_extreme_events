
def periodic(list) :
#Checks for periodicity of extreme events, is an approximation to the method used in the paper.
    if len(list) <= 3 :
        return False
    diff = np.diff(np.array(list))
    median_diff = np.median(diff)
    if median_diff == 0 :
        return False
    ratio = diff / median_diff
    return  np.all((ratio < 1.05) & (ratio > 0.95))

def lower_band(peaks , threshold, extreme_event) :
# Custom function, which finds values which are close to the threshold. 
# If there are more than 3 values in the range of 0.1 of the threshold for 1 extreme event, that event is not considered as an extreme event. 
    almost_extreme_events = [i for i in peaks if i > threshold - 0.1 and i < threshold]
    if len(almost_extreme_events) == 0 :
        return False
    elif len(extreme_event) / len(almost_extreme_events) < 0.3333 :
        return True
    else :
        return False

def num_extreme_events(actual_x_mean) :
#Finds the number of extreme events.
    peaks_indices, _ = find_peaks(actual_x_mean)
    peaks = actual_x_mean[peaks_indices]
    if len(peaks) == 0 :
        return 0
    threshold = max(0.6 , np.mean(peaks) + 8 * np.std(peaks)) #Threshold check
    true = np.array(peaks) > threshold
    extreme_event = peaks[true]
    extreme_event_times = peaks_indices[true]
    if periodic(extreme_event_times) : #Periodicity check
        extreme_event = np.array([])
    if lower_band(peaks , threshold , extreme_event) : #Almost extreme event check
        extreme_event = np.array([])
    if len(extreme_event) / len(peaks) > 0.05 : #Custom critieria, which says that if extreme events are more than 5 percent of the total data then they are not extreme.
        extreme_event = np.array([])
      return len(extreme_event)

  
