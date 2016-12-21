package CoffeeMaker;

public class CoffeeMaker {
    /* Coffee Maker Low Level Hardware Interface */
    enum WarmerPlateStatus { potNotEmpty, potEmpty, warmerEmpty };
    enum WarmerPlateStatus GetWarmerPlateStatus ();
    enum BoilerStatus { boilerEmpty, boilerNotEmpty };
    enum BoilerStatus GetBoilerStatus ();
    enum BrewButtonStatus { brewButtonPushed, brewButtonNotPushed };
    enum BrewButtonStatus GetBrewButtonStatus ();
    enum BoilerState { boilerOn, boilerOff };
    void SetBoilerState (enum BoilerHeaterState s);
    enum WarmerState { warmerOn, warmerOff };
    void SetWarmerState (enum WarmerState s);
    enum IndicatorState { indicatorOn, indicatorOff };
    void SetIndicatorState (enum IndicatorState s);

}
