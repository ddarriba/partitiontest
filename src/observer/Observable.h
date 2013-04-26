/*
 * Observable.hpp
 *
 *  Created on: 30/05/2012
 *      Author: diego
 */

#ifndef OBSERVABLE_H_
#define OBSERVABLE_H_

#include <vector>
#include <string>
#include "Observer.h"
#include "../model/Model.h"

namespace partest {

class Observable {
  public:
    Observable();
    virtual ~Observable();
    void attach(Observer * myObserver);
    void detach(Observer * myObserver);
    void notify_observers(ObservableInfo * info, string * message = 0);
    void notify_observers(MessageType type, unsigned int p_index, Model * model, time_t time = 0,
        unsigned int current_index = 0, unsigned int max_index = 0, string * message = 0);
    void notify_observers(MessageType type, unsigned int p_index, time_t time = 0,
        unsigned int current_index = 0, unsigned int max_index = 0, string * message = 0);
    void notify_observers(MessageType type, unsigned int p_index, ModelSet * modelset, time_t time = 0,
            unsigned int current_index = 0, unsigned int max_index = 0, string * message = 0);
//    void notify_modelset_init(string name, int current_index, int max_index);
  protected:
    void chainObservers(Observable & item);
  private:
    vector<Observer *> observers;
};

} /* namespace partest */
#endif /* OBSERVABLE_H_ */
