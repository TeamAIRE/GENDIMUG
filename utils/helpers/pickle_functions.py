import pickle


def load_pickle(path):
    """
    Deserializes a pickle file.
    :param path: path to the pickle file
    :return: deserialized object
    """
    with open(path, 'rb') as handle:
        return pickle.load(handle)


def save_pickle(obj, path):
    """
    Saves a pickle file.
    :param obj: Object to save
    :param path: str, path to the pickle file
    :return: None
    """
    with open(path, "wb") as f:
        pickle.dump(obj, f)
